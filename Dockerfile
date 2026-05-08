FROM condaforge/miniforge3

RUN apt-get update && apt-get install build-essential git libz-dev -y

RUN git clone https://github.com/xzhuo/mlv.git
WORKDIR /mlv
RUN conda env create --file envs/kmer.yml
RUN conda env create --file envs/align.yml
