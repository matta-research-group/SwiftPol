# Docker file we are building on top of
FROM ghcr.io/mamba-org/micromamba:2.0.8-alpine3.20

# Metadata for the image
LABEL org.opencontainers.image.source="https://github.com/matta-research-group/SwiftPol"
LABEL org.opencontainers.image.description="Tools for building polydisperse polymer systems for molecular dynamics."
LABEL org.opencontainers.image.licenses="BSD-3"

# Set working directory where files will live inside the container
WORKDIR /opt/app

# Copy everything from the current directory (including subdirectories) to /opt/app inside the container
COPY . .

# Fix permissions for the /opt/app directory and subdirectories as root user before switching to mambauser
USER root
RUN chown -R $MAMBA_USER:$MAMBA_USER /opt/app && \
    chmod -R u+rw /opt/app

# Install the environment and Jupyter (meta-package includes JupyterLab and Jupyter Notebook)
RUN micromamba install -y -n base -f ./Dev_tools/docker.yml && \
    micromamba install -y -n base -c conda-forge jupyter 
    #&& \
    #micromamba clean --all --yes


# Enable environment for future RUN/CMDs
ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN pip install ./

# Ensure mount directory exists (useful for volumes)
RUN mkdir -p /opt/app
USER $MAMBA_USER

# Expose port 8888 (for Jupyter)
EXPOSE 8888

# Optional: Set the default command to start Jupyter Notebook
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--allow-root", "--no-browser"]



