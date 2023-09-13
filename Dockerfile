

FROM python:3.8-slim-buster

ADD scripts/get_FDA_thresholds.py /opt/scripts/get_FDA_thresholds.py
ADD scripts/get_neoantigen_qc.py /opt/scripts/get_neoantigen_qc.py
ADD scripts/requirements.txt /opt/scripts/requirements.txt
ADD scripts/fda_quality_thresholds.csv /opt/scripts/fda_quality_thresholds.csv
ADD scripts/generate_reviews_files.py /opt/scripts/generate_reviews_files.py
ADD scripts/bold_classII.py /opt/scripts/bold_classII.py
ADD scripts/color_peptides51mer.py /opt/scripts/color_peptides51mer.py
ADD scripts/setup_reivew.py /opt/scripts/setup_review.py

RUN chmod +r /opt/scripts/*


RUN pip install --upgrade pip


RUN pip3 install -r /opt/scripts/requirements.txt
