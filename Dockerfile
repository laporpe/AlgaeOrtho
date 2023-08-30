FROM continuumio/anaconda3:2022.05

# set the working directory
WORKDIR /app

# copy the environment file first, to avoid re-running conda install on every code change
COPY environment.yml ./
RUN conda env create -f environment.yml

# install unzip
RUN apt-get update && apt-get install -y unzip

# make RUN commands use the new environment
# SHELL ["conda", "run", "-n", "algaeortho", "/bin/bash", "-c"]
ENV PATH="/opt/conda/envs/algaeortho/bin:$PATH"
RUN which python && python --version

# copy the zip file and unzip it
COPY ortholog_groups.tsv.zip ./
RUN unzip ortholog_groups.tsv.zip

# copy the rest of the files
COPY . .

# run the app without starting the server
# to generate the df_to_merge pickle file
RUN python -u main.py --noserver

# set environment variables
ENV SERVER_PORT=8050
ENV DASH_DEBUG_MODE=True

# expose the server port
EXPOSE 8050

### everything above this runs only when building the image
### everything below this runs when the container is started

# run the app
CMD ["python", "-u", "main.py"]