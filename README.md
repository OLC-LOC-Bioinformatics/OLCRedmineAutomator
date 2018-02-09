# OLCRedmineAutomator

The OLCRedmineAutomator allows for easy access to a number of
bioinformatics tools for use through Redmine.

The majority of our tools work by parsing a list of Sample IDs provided
in the description field of a Redmine issue submitted through the
[OLC Redmine portal](http://redmine.biodiversity.agr.gc.ca/projects/cfia/).
The requested tool to run the analysis must be specified in the subject line.

Below are a list of acceptable keywords that the OLCRedmineAutomator
will detect in an issue subject line. The amount of resources allocated
to the OLC Slurm cluster for each respective job type is also displayed.

| Keyword          | CPUs |  RAM (GB)|
| ---------------  |:----:|:--------:|
| Strainmash       | 8    |  12      |
| WGS Assembly     | 12   |  192     |
| AutoCLARK        | 48   |  192     |
| SNVPhyl          | 48   |  192     |
| Diversitree      | 56   |  192     |
| External Retrieve| 1    |  6       |
| PlasmidExtractor | 48   |  192     |
| Merge            | 48   |  192     |
| CloseRelatives   | 48   |  192     |


### Tool Descriptions (Work in progress)
Unless otherwise stated, the **Usage** section of each tool description
provides an example of the input to provide in a newly created Redmine
issue's description field. The name of the requested tool must be
specified in the subject line of the issue.

#### Strainmash
**Description:** Reads a list of SeqIDs from a Redmine issue and
calls Mash for each against a sketch of the entire GenBank
strain database (Up to date as of Dec. 01, 2017).
Returns a formatted Mash.screen output file per SeqID.

**Usage:**
```
2017-SEQ-0918
2017-SEQ-0919
2017-SEQ-0920
2017-SEQ-0921
```

---
#### Diversitree
**Description:** Takes an integer *n* and a list of SeqIDs from a Redmine issue.
This program will select the *n* most representative strains from the given list of SeqIDs.

**Usage:**
```
2
2017-SEQ-0918
2017-SEQ-0919
2017-SEQ-0920
2017-SEQ-0921
```

---
#### AutoCLARK
**Description:**

**Usage:**

---
#### CloseRelatives
**Description:**

**Usage:**

---
#### PlasmidExtractor
**Description:** Takes a list of SeqIDs and runs PlasmidExtractor on each sample.
For details on the method, see __https://lowandrew.github.io/Plasmid_Assembler/__

**Usage:**
```
2017-SEQ-0918
2017-SEQ-0919
2017-SEQ-0920
2017-SEQ-0921
```

---
#### WGSAssembly
**Description:**

**Usage:**

---
#### SNVPHyl
**Description:**

**Usage:**

---


### Internal Notes
1. Log into the head node (ubuntu@192.168.1.26)

2. Activate the virtual environment
    - ```source /mnt/nas/Redmine/.virtualenvs/OLCRedmineAutomator/bin/activate```

3. Call automation script
    - ```python /mnt/nas/Redmine/OLCRedmineAutomator/api.py 2> /dev/null```

4. Enjoy Redmine automation

Running this program requires a *settings.py* file that is not included in
this repository for security. Here's a censored example of setup.py:

```
API_KEY = 'super_secret'
BIO_REQUESTS_DIR = '/mnt/nas/bio_requests'
AUTOMATOR_KEYWORDS = {
    'strainmash':
        {'memory': 12000, 'n_cpu': 8},
    'autoclark':
        {'memory': 192000, 'n_cpu': 48},
    'diversitree':
        {'memory': 192000, 'n_cpu': 56},
}
```