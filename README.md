# AlgaeOrtho Dash App
### Author: Mary-Francis LaPorte, Neha Arora, Struan Clark, Ambarish Nag


## Description
This Dash app, named AlgaeOrtho, is designed to analyze and visualize orthologous gene groups in algae. It provides functionality to load data, perform Clustal Omega alignment, and generate phylogenetic trees.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Step-by-step Guide](#guide)

## Installation
To run this Dash app using Docker without cloning and building the repository locally, run:
- Mac/Linux/WSL:
  - `docker run -p 8050:8050 -v ~/.algaeortho_data:/data ghcr.io/laporpe/algaeortho:main`
- Windows:
  - `docker run -p 8050:8050 -v "%HOMEPATH%\.algaeortho_data":/data ghcr.io/laporpe/algaeortho:main`

## Installation (from source)
To build and run this Dash app using Docker Compose, follow these steps:
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

## Guide

### Helpful software and pre-requisites
1.	Docker: a tool to handle the software prerequisites for the AlgaeOrtho application: https://hub.docker.com/. Once Docker is installed, it will handle the Python software installation that our application uses to run. Although it will require the user to run a couple of lines of code in the terminal, this is easier than the user needing to handle software installation and versioning themselves. 

2.	GitHub: a tool to share code and make it accessible over the internet. The user will want a GitHub account https://github.com/  to retrieve the code for this application: https://github.com/laporpe/AlgaeOrtho . This code can then be downloaded to the user’s device so they can run the code. 

3.	One way to download the code for the application is through the GitHub Desktop application: https://desktop.github.com/ . The GitHub Desktop application can be connected with the users’ GitHub account, so that the user can run AlgaeOrtho code (or any other code accessed from GitHub) on their computer. Here are instructions on how to do that: https://docs.github.com/en/desktop/adding-and-cloning-repositories/cloning-a-repository-from-github-to-github-desktop 
![desktopapp](https://github.com/laporpe/AlgaeOrtho/assets/17015641/0bec2337-36f1-44c6-9746-95cccd79518d)
4.	Visual Studio Code (VS Code): a code editor that allows the user to open up the AlgaeOrtho files and run the application. VS Code is advantageous because it has a built-in terminal. Once the code has been obtained from Github, the user should open it in VS Code.
![vscode](https://github.com/laporpe/AlgaeOrtho/assets/17015641/2ef774c9-860a-4aae-aded-c38fe47aa71b)

### How to use Docker?

1.	Download Docker https://hub.docker.com/ 
2.	Pull the code for the application from Github https://github.com/laporpe/AlgaeOrtho. 
a.	 One way to do this is through using the GitHub Desktop App: https://desktop.github.com/. See link in Part 3 above. 
 
3.	Open the directory that contains the code for the application locally in VS Code. This can be done by pressing the “Open in Visual Studio Code” button in the GitHub Desktop application (center section, two boxes beneath the “No local changes” heading). 
4.	Once the application files have been opened in VS Code, navigate to the Terminal in VS code (usually at the bottom of the screen), and  run this: 

docker compose build 

a.	This command builds the docker container so that all of the necessary prerequisites for this application are pre-installed. 
 
5.	Once that command has finished, In the Terminal in VS code (with the application directory open), run this command: 
docker compose up
![openapp](https://github.com/laporpe/AlgaeOrtho/assets/17015641/aa9b3e9f-b948-4fc1-bf1f-20852924ad6f)
 
6.	This will make the application accessible locally at http://127.0.0.1:8050 (another equivalent option: localhost:8050)
a.	In other words, type 127.0.0.1:8050 (or localhost:8050) into your browser (for example, Google Chrome) and the application will open. The application is running locally, and using the graphical tools from your browser to display the application. As the application is running locally, anything you upload only is being processed on your computer and is not being uploaded online. 
7.	From here you should be able to use the example buttons as well as upload your own sequences!

### How to use the Algae Ortho application?
 
This is how the application will appear when it is first loaded. Select the buttons on the upper left to visualize pre-loaded data (described in the manuscript). 

![app](https://github.com/laporpe/AlgaeOrtho/assets/17015641/1b586f85-e64a-4b67-9a14-d3178287e943)

To upload your own data, use the “upload your protein query” button. Upload any protein .fasta file, that does not end each protein with asterisks (*). Try a short file with only a few proteins first, as these run quickly. 
Leave the checkbox “unchecked” for the “quick” version of the application, which only finds the putative orthologs from the sonic paranoid results, but does not visualize them.  This will allow you to download the resulting protein fasta file using the “download results” button. Leave the checkbox “checked” to visualize the results of the putative ortholog file in the application. 
Download all results files by clicking the “download results” button. This will download a folder that needs to be unzipped, that will contain the.fasta file of all putative orthologs, as well as the generated Newick file to visualize a tree of those sequences in your tree visualizer of choice. 
Remember to scroll down (on the white-background section of the application) to see the ortholog tree!  
The figures can be zoomed in and out using the Dash app controls  
<img width="262" alt="dashappcontrols" src="https://github.com/laporpe/AlgaeOrtho/assets/17015641/b60f4b49-56ca-4c72-b5c8-306df200ab4b">

### Troubleshooting: 
If you are running this on an M1 Mac: the visualization aspect of the application uses a function that cannot work with an M1 Mac. 

There is a separate container distributed for users on an M1 Mac, that has a “checkbox unchecked” version of the application. To access this container, follow the instructions above, but replace the docker commands with the following: 

docker compose -f docker-compose-noclustalo.yml build
docker compose -f docker-compose-noclustalo.yml up


This specifies the “noclustalo” (non-visualized) version, which should work on all machines.

