FROM python:2.7.15-slim

COPY requirements.txt /tmp/
RUN cd /tmp/ && pip install -r requirements.txt && rm /tmp/requirements.txt

COPY *.py /work/
WORKDIR /work

CMD ["/bin/bash"]
