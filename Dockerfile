FROM python:3.7.3-stretch
WORKDIR /varlap
COPY . .

# Install the python package (and executable)
RUN pip3 install .
