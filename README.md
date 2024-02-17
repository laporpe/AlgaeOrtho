# AlgaeOrtho Dash App
### Author: M. F. Laporte

## Description
This Dash app, named AlgaeOrtho, is designed to analyze and visualize orthologous gene groups in algae. It provides functionality to load data, perform Clustal Omega alignment, and generate phylogenetic trees.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)

## Installation
To run this Dash app using Docker Compose, follow these steps:
1. Clone this repository.
2. Navigate to the project directory.
3. Build and start the Docker container using the command: `docker-compose up --build`
4. Once the containers are up and running, access the app through your web browser at `http://localhost:8050/`.

## Usage
Once the app is running, you can access it through your web browser at `http://localhost:8050/`. The app allows you to perform the following tasks:

### Select Input
- Choose pre-loaded data by clicking on the provided buttons.
- Alternatively, upload your own protein query file.

### Configure Clustal Omega
- Enable/disable the option to run Clustal Omega for alignment.

### Download Results
- Download the generated results, including alignment files and percent identity matrices (PIM).
- Several result types are available based on the context of the analysis.

#### Distance Matrix Files:
- **Mode**: These files are provided when the application is run with the pre-loaded data.
- **Description:** These files contain the percent identity matrix (PIM) calculated from the aligned sequences.
- **Format:** Text files (`.pim.txt`).

#### Merged Fasta Files:
- **Mode**: These files are provided when the application is run with user data and Clustal Omega disabled.
- **Description:** These files contain merged fasta sequences from user-uploaded data and pre-loaded datasets.
- **Format:** Fasta format (`.fasta` files).

#### Results Zip Files:
- **Mode**: These files are provided when the application is run with user data and Clustal Omega enabled.
- **Description:** These files contain a compressed archive of all relevant output files, including alignment files, distance matrix files, and merged fasta files.
- **Format:** ZIP archive (`.zip` files).

### Additional Docker Commands

- Build: `docker compose build`
- Run: `docker compose up`
- Cleanup: `docker compose down --remove-orphans`
