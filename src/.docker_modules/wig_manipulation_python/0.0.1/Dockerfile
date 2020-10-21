FROM python:3.7-alpine

COPY wig2bed.py /usr/bin/
COPY wig_average_normalization.py /usr/bin/
RUN apk add  --update --no-cache bash && \
  chmod +x /usr/bin/wig2bed.py /usr/bin/wig_average_normalization.py
