# I can't get it to work :(

FROM python:3.8-slim-buster

ADD scripts/get_FDA_thresholds.py /opt/scripts/get_FDA_thresholds.py
ADD scripts/get_neoantigen_qc.py /opt/scripts/get_neoantigen_qc.py
ADD scripts/requirements.txt /opt/scripts/requirements.txt
ADD scripts/fda_quality_thresholds.csv /opt/scripts/fda_quality_thresholds.csv


RUN chmod +r /opt/scripts/get_FDA_thresholds.py
RUN chmod +r /opt/scripts/get_neoantigen_qc.py
RUN chmod +r /opt/scripts/fda_quality_thresholds.csv

RUN pip install --upgrade pip


RUN pip3 install -r /opt/scripts/requirements.txt
