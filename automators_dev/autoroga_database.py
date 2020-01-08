import sqlalchemy as sa
import datetime
import time

from automator_settings import POSTGRES_PASSWORD, POSTGRES_USERNAME


# Default connection to the address of the head node - should be 192.168.1.5
def connect(user, password, db, host='192.168.1.5', port=5432):
    url = 'postgresql://{}:{}@{}:{}/{}'
    url = url.format(user, password, host, port, db)

    # To prevent race conditions, the POSTGRES_USERNAME should be set up with a connection limit of 1.
    # Then, if other users are currently accessing ROGA database, just wait.
    connection_created = False

    while connection_created is False:
        try:
            # The return value of create_engine() is the connection object
            con = sa.create_engine(url, client_encoding='utf8')

            # Bind the connection to MetaData()
            meta = sa.MetaData(bind=con, reflect=True)
            connection_created = True
        except sa.exc.OperationalError:  # This is the error you get if too many users try to connect.
            time.sleep(10)

    return con, meta


def update_db(date, year, genus, lab, source, amendment_flag, amended_id):
    con, meta = connect(user=POSTGRES_USERNAME, password=POSTGRES_PASSWORD, db='autorogadev')  # NOTE db=autorogadev

    ROGA_ID_SEQ = sa.Sequence('roga_id_seq')

    try:  # Create table if it doesn't already exist
        autoroga_project_table = sa.Table('autoroga_project_table', meta,
                                          sa.Column('id', sa.INTEGER, ROGA_ID_SEQ,
                                                    primary_key=True, server_default=ROGA_ID_SEQ.next_value()),
                                          sa.Column('roga_id', sa.String(64)),
                                          sa.Column('genus', sa.String(64)),
                                          sa.Column('lab', sa.String(16)),
                                          sa.Column('source', sa.String(64)),
                                          sa.Column('amendment_flag', sa.String(16)),
                                          sa.Column('amended_id', sa.String(64)),
                                          sa.Column('date', sa.Date),
                                          sa.Column('time', sa.DateTime, default=datetime.datetime.utcnow),
                                          sa.Column('deletion_date', sa.Date),
                                          sa.Column('deletion_reason', sa.String(256))
                                          )
        meta.create_all()
        print('Successfully created autoroga_project_table')

    except:  # Retrieve table if it already exists
        autoroga_project_table = sa.Table('autoroga_project_table', meta, autoload=True, autoload_with=sa.engine)
        print('Successfully retrieved autoroga_project_table')

    # Grab what the next key value will be
    select_next_value = sa.select([autoroga_project_table.c.id])
    keys = con.execute(select_next_value)

    try:
        next_val = max(keys)[0] + 1
    except:
        next_val = 1

    # Create ROGA ID
    select_next_roga_id = sa.select([autoroga_project_table.c.roga_id])
    keys = con.execute(select_next_roga_id)
    roga_ids = keys.fetchall()
    # Now parse through ROGA IDs to figure out what the next one should be. roga_ids are a list of tuples
    ids_for_year = list()
    for item in roga_ids:
        # ROGA is actually first element of tuple.
        roga = item[0]
        roga_year = int(roga.split('-')[0])
        roga_id = int(roga.split('-')[-1])
        if roga_year == int(year):
            ids_for_year.append(roga_id)

    i = 1
    roga_id_found = False
    while roga_id_found is False:
        if i in ids_for_year:
            i += 1
        else:
            roga_id_found = True

    roga_id = str(year) + '-ROGA-DEV-' + '{:04d}'.format(i)

    # Insert new row into autoroga_project_table table
    ins = autoroga_project_table.insert().values(roga_id=roga_id, genus=genus, date=date, lab=lab, source=source,
                                                 amendment_flag=amendment_flag, amended_id=amended_id,
                                                 time=datetime.datetime.utcnow())
    con.execute(ins)

    return roga_id
