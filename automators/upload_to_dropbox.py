"""
This script uploads a specified local file to Dropbox using the Dropbox API.

Usage:
    python upload_to_dropbox.py -f /path/to/local/file -d /path/in/dropbox/file
    python upload_to_dropbox.py -c AUTHORIZATION_CODE

Arguments:
    -f, --file: Path to the local file to be uploaded.
    -d, --destination: Full path in Dropbox where the file will be uploaded,
                       including the file name.
    -c, --code: Authorization code to obtain initial access and refresh tokens.

The Dropbox access token and refresh token are imported from the tokens.py file
"""
# Standard imports
import argparse
import os
import time

# Third-party imports
import dropbox
from dropbox.dropbox_client import Dropbox
import requests
from tqdm import tqdm

# Local imports
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_REFRESH_TOKEN,
    DROPBOX_APP_KEY,
    DROPBOX_APP_SECRET
)

# Define constants for chunk size, timeout, and maximum retries
CHUNK_SIZE = 50 * 1024 * 1024  # 50MB chunks
TIMEOUT = 1500  # 25 minutes
MAX_RETRIES = 5  # Maximum number of retries for each chunk

class DropboxWithTimeout(Dropbox):
    """
    Subclass of Dropbox client with a custom timeout setting.
    """
    def __init__(
        self, 
        oauth2_access_token, 
        timeout=TIMEOUT, 
        **kwargs
    ):
        super().__init__(oauth2_access_token, **kwargs)
        self._session.timeout = timeout

def refresh_access_token(
    refresh_token: str, 
    app_key: str, 
    app_secret: str
) -> str:
    """
    Refresh the Dropbox access token using the refresh token.

    :param refresh_token: str, Dropbox API refresh token
    :param app_key: str, Dropbox API app key
    :param app_secret: str, Dropbox API app secret
    :return: str, New access token
    """
    url = "https://api.dropbox.com/oauth2/token"
    data = {
        "grant_type": "refresh_token",
        "refresh_token": refresh_token,
        "client_id": app_key,
        "client_secret": app_secret,
    }
    response = requests.post(url, data=data)
    response_data = response.json()
    if "access_token" in response_data:
        return response_data["access_token"]
    else:
        raise Exception("Failed to refresh access token")

def get_initial_tokens(
    app_key: str, 
    app_secret: str, 
    authorization_code: str
) -> tuple:
    """
    Obtain initial access and refresh tokens using the authorization code.

    :param app_key: str, Dropbox API app key
    :param app_secret: str, Dropbox API app secret
    :param authorization_code: str, Authorization code obtained from Dropbox
    :return: tuple, Access token and refresh token
    """
    url = "https://api.dropbox.com/oauth2/token"
    data = {
        "code": authorization_code,
        "grant_type": "authorization_code",
        "client_id": app_key,
        "client_secret": app_secret,
    }
    response = requests.post(url, data=data)
    response_data = response.json()
    if "access_token" in response_data and "refresh_token" in response_data:
        return response_data["access_token"], response_data["refresh_token"]
    else:
        print(response_data)
        raise Exception("Failed to obtain initial tokens")
    

def update_settings_file(
    access_token: str, 
    refresh_token: str
) -> None:
    """
    Update the tokens.py file with new access and refresh tokens.

    :param access_token: str, New access token
    :param refresh_token: str, New refresh token
    """
    # Set the path of the tokens.py file
    tokens_file_path = os.path.join(
        '/mnt/nas2/redmine/applications/OLCRedmineAutomator/automators', 
        'tokens.py'
    )
    
    with open(tokens_file_path, 'w', encoding='utf-8') as token_file:
        token_file.write('DROPBOX_ACCESS_TOKEN = "{access_token}"\n'.format(
            access_token=access_token)
        )
        token_file.write('DROPBOX_REFRESH_TOKEN = "{refresh_token}"\n'.format(
            refresh_token=refresh_token)
        )
        token_file.write('DROPBOX_APP_KEY = "{DROPBOX_APP_KEY}"\n'.format(
            DROPBOX_APP_KEY=DROPBOX_APP_KEY)
        )
        token_file.write(
            'DROPBOX_APP_SECRET = "{DROPBOX_APP_SECRET}"\n'.format(
                DROPBOX_APP_SECRET=DROPBOX_APP_SECRET
            )
        )


def check_and_refresh_access_token(
    access_token: str,
    refresh_token: str,
    app_key: str,
    app_secret: str
) -> str:
    """
    Check if the access token is expired and refresh it if necessary.

    :param access_token: str, Dropbox API access token
    :param refresh_token: str, Dropbox API refresh token
    :param app_key: str, Dropbox API app key
    :param app_secret: str, Dropbox API app secret
    :return: str, Valid access token
    """
    dbx = DropboxWithTimeout(access_token)
    try:
        # Try to get account info to check if the access token is valid
        dbx.users_get_current_account()
        return access_token
    except dropbox.exceptions.AuthError as exc:
        if exc.error.is_expired_access_token():
            print("Access token expired. Refreshing access token...")
            new_access_token = refresh_access_token(
                refresh_token, app_key, app_secret
            )
            update_settings_file(new_access_token, refresh_token)
            return new_access_token
        else:
            raise


def upload_to_dropbox(
    access_token: str, 
    refresh_token: str, 
    app_key: str,
    app_secret: str, 
    local_file_path: str,
    dropbox_destination_path: str = '/Carling-FTP'
) -> str:
    """
    Uploads a file to Dropbox and generates a shared link.

    :param access_token: str, Dropbox API access token
    :param refresh_token: str, Dropbox API refresh token
    :param app_key: str, Dropbox API app key
    :param app_secret: str, Dropbox API app secret
    :param local_file_path: str, Path to the local file to be uploaded
    :param dropbox_destination_path: str, Full path in Dropbox where the file
        will be uploaded, including the file name
    :return: str, Shared link to the uploaded file
    """
    # Ensure that the local file exists
    if not check_local_file_exists(local_file_path):
        return ""

    # Ensure that the Dropbox destination path is valid
    dropbox_destination_path = prepare_dropbox_destination_path(
        dropbox_destination_path=dropbox_destination_path,
        local_file_path=local_file_path
    )
    
    # Check if the access token is expired and refresh it if necessary
    access_token = check_and_refresh_access_token(
        access_token=access_token,
        refresh_token=refresh_token,
        app_key=app_key,
        app_secret=app_secret
    )
    
    # Create a Dropbox client
    dbx = DropboxWithTimeout(access_token)
    
    # Clear out files that were uploaded more than 14 days ago
    clear_old_files(dbx=dbx)

    # Check to see if the file is already in Dropbox
    if check_file_exists_in_dropbox(
        dbx=dbx,
        dropbox_destination_path=dropbox_destination_path
    ):
        return ""

    # Upload the file to Dropbox
    upload_file(
        dbx=dbx,
        local_file_path=local_file_path,
        dropbox_destination_path=dropbox_destination_path
    )
    print(
        "File {local_file_path} uploaded to Dropbox at "
        "{dropbox_destination_path}".format(
            local_file_path=local_file_path,
            dropbox_destination_path=dropbox_destination_path
        )
    )
    
    # Create a shared link for the uploaded file
    return create_shared_link(
        dbx=dbx,
        dropbox_destination_path=dropbox_destination_path
    )


def clear_old_files(dbx: Dropbox) -> None:
    """
    Clear out any files in the Dropbox account that are older than 14 days.

    :param dbx: Dropbox client
    """
    
    # Calculate the timestamp for 14 days ago
    fourteen_days_ago = time.time() - 14 * 24 * 60 * 60
    try:
        result = dbx.files_list_folder('', recursive=True)
        for entry in result.entries:

            # Ensure that the entry contains a FileMetadata object
            if isinstance(entry, dropbox.files.FileMetadata):

                # Check if the file was last modified more than 14 days ago
                if entry.server_modified.timestamp() < fourteen_days_ago:
                    dbx.files_delete_v2(entry.path_lower)
                    print(
                        "Deleted old file: {path}".format(
                            path=entry.path_lower
                        )
                    )
    except dropbox.exceptions.ApiError as exc:
        print("Error clearing old files: {exc}".format(exc=exc))


def check_local_file_exists(local_file_path: str) -> bool:
    """
    Check if the local file exists.

    :param local_file_path: str, Path to the local file
    :return: bool, True if the file exists, False otherwise
    """
    if not os.path.isfile(local_file_path):
        print("Error: The file {local_file_path} does not exist.".format(
            local_file_path=local_file_path
        ))
        return False
    return True

def prepare_dropbox_destination_path(
    dropbox_destination_path: str, local_file_path: str
) -> str:
    """
    Prepare the Dropbox destination path.

    :param dropbox_destination_path: str, Dropbox destination path
    :param local_file_path: str, Local file path
    :return: str, Prepared Dropbox destination path
    """
    if not dropbox_destination_path.startswith('/'):
        dropbox_destination_path = '/' + dropbox_destination_path

    dropbox_destination_path = dropbox_destination_path.rstrip('/') + '/'
    dropbox_destination_path += os.path.basename(local_file_path)
    return dropbox_destination_path

def check_file_exists_in_dropbox(
    dbx: Dropbox, dropbox_destination_path: str
) -> bool:
    """
    Check if the file already exists in Dropbox.

    :param dbx: Dropbox client
    :param dropbox_destination_path: str, Dropbox destination path
    :return: bool, True if the file exists, False otherwise
    """
    try:
        dbx.files_get_metadata(dropbox_destination_path)
        print(
            "Error: The file already exists at {dropbox_path}".format(
                dropbox_path=dropbox_destination_path
            )
        )
        return True
    except dropbox.exceptions.ApiError as exc:
        if exc.error.is_path() and exc.error.get_path().is_not_found():
            return False
        else:
            raise

def upload_file_in_chunks(
    dbx: Dropbox, f, dropbox_destination_path: str, file_size: int
) -> None:
    """
    Upload a file to Dropbox in chunks.

    :param dbx: Dropbox client
    :param f: File object
    :param dropbox_destination_path: str, Dropbox destination path
    :param file_size: int, Size of the file
    """
    upload_session_start_result = dbx.files_upload_session_start(
        f.read(CHUNK_SIZE)
    )
    cursor = dropbox.files.UploadSessionCursor(
        session_id=upload_session_start_result.session_id, offset=f.tell()
    )
    commit = dropbox.files.CommitInfo(path=dropbox_destination_path)

    progress_bar = tqdm(
        total=file_size, unit='B', unit_scale=True, desc='Uploading', ncols=80
    )
    while f.tell() < file_size:
        if (file_size - f.tell()) <= CHUNK_SIZE:
            retry_count = 0
            while retry_count < MAX_RETRIES:
                try:
                    dbx.files_upload_session_finish(
                        f.read(CHUNK_SIZE), cursor, commit
                    )
                    progress_bar.update(CHUNK_SIZE)
                    break
                except dropbox.exceptions.AuthError:
                    access_token = refresh_access_token(
                        refresh_token=DROPBOX_REFRESH_TOKEN,
                        app_key=DROPBOX_APP_KEY,
                        app_secret=DROPBOX_APP_SECRET
                    )
                    dbx = DropboxWithTimeout(access_token)
                except Exception as exc:
                    print("Error: {exc}. Retrying...".format(exc=exc))
                    retry_count += 1
                    time.sleep(2 ** retry_count)
        else:
            retry_count = 0
            while retry_count < MAX_RETRIES:
                try:
                    dbx.files_upload_session_append_v2(
                        f.read(CHUNK_SIZE), cursor
                    )
                    cursor.offset = f.tell()
                    progress_bar.update(CHUNK_SIZE)
                    break
                except dropbox.exceptions.AuthError:
                    access_token = refresh_access_token(
                        refresh_token=DROPBOX_REFRESH_TOKEN,
                        app_key=DROPBOX_APP_KEY,
                        app_secret=DROPBOX_APP_SECRET
                    )
                    dbx = DropboxWithTimeout(access_token)
                except Exception as exc:
                    print("Error: {exc}. Retrying...".format(exc=exc))
                    retry_count += 1
                    time.sleep(2 ** retry_count)
    progress_bar.close()


def upload_file(
    dbx: Dropbox, local_file_path: str, dropbox_destination_path: str
) -> None:
    """
    Upload a file to Dropbox.

    :param dbx: Dropbox client
    :param local_file_path: str, Local file path
    :param dropbox_destination_path: str, Dropbox destination path
    """
    with open(local_file_path, 'rb') as f:
        file_size = os.path.getsize(local_file_path)
        if file_size <= CHUNK_SIZE:
            dbx.files_upload(f.read(), dropbox_destination_path)
        else:
            upload_file_in_chunks(dbx, f, dropbox_destination_path, file_size)


def create_shared_link(
    dbx: Dropbox, dropbox_destination_path: str
) -> str:
    """
    Create a shared link for the uploaded file.

    :param dbx: Dropbox client
    :param dropbox_destination_path: str, Dropbox destination path
    :return: str, Direct download link
    """
    try:
        shared_link_metadata = dbx.sharing_create_shared_link_with_settings(
            dropbox_destination_path
        )
        shared_link = shared_link_metadata.url
        direct_download_link = shared_link.replace(
            "www.dropbox.com", "dl.dropboxusercontent.com"
        ).replace("dl=0", "dl=1")
        print("Shared link: {shared_link}".format(shared_link=shared_link))
        return direct_download_link
    except dropbox.exceptions.ApiError as exc:
        print("Error creating shared link: {exc}".format(exc=exc))
        return ""



def main() -> None:
    """
    Main function to parse arguments and upload the file to Dropbox.
    """
    parser = argparse.ArgumentParser(description='Upload a file to Dropbox.')
    parser.add_argument(
        '-f', '--file',
        type=str,
        help='Path to the local file to be uploaded'
    )
    parser.add_argument(
        '-d', '--destination',
        type=str,
        required=False,
        default='/Carling-FTP',
        help='Full path in Dropbox where the file will be uploaded. '
        'Must start with a "/", but a "/" will be added if you don\'t provide '
        'one. Default value is "/Carling-FTP".'
    )
    parser.add_argument(
        '-c', '--code',
        type=str,
        required=False,
        help=(
            'Authorization code to obtain initial access and refresh tokens. '
            'To get this code, follow these steps:\n'
            '1. Generate the authorization URL:\n'
            '   https://www.dropbox.com/oauth2/authorize?client_id='
            '<APP_KEY>&response_type=code&redirect_uri=<REDIRECT_URI>\n'
            '2. Open the URL in a web browser and authorize the app.\n'
            '3. Extract the authorization code from the redirect URL.'
        )
    )
    parser.add_argument(
        '-g', '--generate-auth-url',
        required=False,
        action='store_true',
        help='Generate the authorization URL to obtain the authorization code'
    )  

    args = parser.parse_args()

    # Generate the authorization URL
    if args.generate_auth_url:
        authorization_url = (
            "https://www.dropbox.com/oauth2/authorize?"
            "client_id={app_key}&response_type=code&token_access_type=offline"
        ).format(app_key=DROPBOX_APP_KEY)
        
        print("Authorization URL: {authorization_url}".format(
            authorization_url=authorization_url)
            )
        raise SystemExit
    if args.code:
        # Obtain initial tokens using the authorization code
        access_token, refresh_token = get_initial_tokens(
            app_key=DROPBOX_APP_KEY,
            app_secret=DROPBOX_APP_SECRET,
            authorization_code=args.code
        )
        # Update the settings file with the new tokens
        update_settings_file(access_token, refresh_token)
        print("Initial tokens obtained and tokens.py updated.")
    elif args.file:
        # Upload the file to Dropbox and get the shared link
        shared_link = upload_to_dropbox(
            access_token=DROPBOX_ACCESS_TOKEN,
            refresh_token=DROPBOX_REFRESH_TOKEN,
            app_key=DROPBOX_APP_KEY,
            app_secret=DROPBOX_APP_SECRET,
            local_file_path=args.file,
            dropbox_destination_path=args.destination
        )
        if shared_link:
            print("Download link: {shared_link}".format(
                shared_link=shared_link)
            )
    else:
        # Print help message if arguments are missing
        parser.print_help()

if __name__ == "__main__":
    main()
