FROM ncbi/blast:2.7.1 as blast

FROM python:3.5.7-slim

COPY --from=blast /blast/ /blast/
ENV BLASTDB /blast/blastdb
ENV PATH="${PATH}:/blast/bin:/usearch/bin"
COPY requirements.txt /tmp/
RUN apt-get update && apt-get install -y gcc libgomp1 && \
	cd /tmp/ && pip install -r requirements.txt && \
	rm /tmp/requirements.txt && \
	apt-get remove --purge -y gcc

COPY deps/ /work/deps/
COPY Snakefile *.py *.sh /work/

WORKDIR /work
CMD ["/bin/bash"]
