FROM python:3.7-alpine

COPY check_design_chip.py /usr/bin/
COPY check_design_chip_quant.py /usr/bin/
RUN apk add  --update --no-cache bash && \
  chmod +x /usr/bin/check_design_chip.py && \
  chmod +x /usr/bin/check_design_chip_quant.py
