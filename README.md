# libb_rspreprocessing

## Tools used in this project
* [hydra](https://hydra.cc/): Manage configuration files - [article](https://mathdatasimplified.com/stop-hard-coding-in-a-data-science-project-use-configuration-files-instead/)
* [pdoc](https://github.com/pdoc3/pdoc): Automatically create an API documentation for your project
* [pre-commit plugins](https://pre-commit.com/): Automate code reviewing formatting


## Project Structure

```bash
.
├── config                      
│   ├── main.yaml                   # Main configuration file
│   ├── model                       # Configurations for training model
│   │   ├── model1.yaml             # First variation of parameters to train model
│   │   └── model2.yaml             # Second variation of parameters to train model
│   └── process                     # Configurations for processing data
│       ├── process1.yaml           # First variation of parameters to process data
│       └── process2.yaml           # Second variation of parameters to process data
├── data            
│   ├── final                       # data after training the model
│   ├── processed                   # data after processing
│   └── raw                         # raw data
├── docs                            # documentation for your project
├── .gitignore                      # ignore files that cannot commit to Git
├── Makefile                        # store useful commands to set up the environment
├── models                          # store models
├── notebooks                       # store notebooks
├── pyproject.toml                  # Configure black

├── README.md                       # describe your project
├── src                             # store source code
│   ├── __init__.py                 # make src a Python module 
│   ├── process.py                  # process data before training model
│   ├── train_model.py              # train model
│   └── utils.py                    # store helper functions
└── tests                           # store tests
    ├── __init__.py                 # make tests a Python module 
    ├── test_process.py             # test functions for process.py
    └── test_train_model.py         # test functions for train_model.py
```

## Set up the environment


1. Create the virtual environment:
```bash
python3 -m venv venv
```
2. Activate the virtual environment:

- For Linux/MacOS:
```bash
source venv/bin/activate
```
- For Command Prompt:
```bash
.\venv\Scripts\activate
```
3. Install dependencies:
- To install all dependencies, run:
```bash
pip install -r requirements-dev.txt
```
- To install only production dependencies, run:
```bash
pip install -r requirements.txt
```
- To install a new package, run:
```bash
pip install <package-name>
```


## View and alter configurations
To view the configurations associated with a Pythons script, run the following command:
```bash
python src/process.py --help
```
Output:
```yaml
process is powered by Hydra.

== Configuration groups ==
Compose your configuration from those groups (group=option)

model: model1, model2
process: process1, process2


== Config ==
Override anything in the config (foo.bar=value)

process:
  use_columns:
  - col1
  - col2
model:
  name: model1
data:
  raw: data/raw/sample.csv
  processed: data/processed/processed.csv
  final: data/final/final.csv
```

To alter the configurations associated with a Python script from the command line, run the following:
```bash
python src/process.py data.raw=sample2.csv
```

## Auto-generate API documentation

To auto-generate API document for your project, run:

```bash
make docs
```
