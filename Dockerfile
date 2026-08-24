# Use a lightweight Python base image
FROM python:3.12-slim

# Set working directory inside the container
WORKDIR /app

# Copy only the server requirements first to leverage Docker cache layers
COPY server/requirements.txt /app/requirements.txt

# Install the server dependencies
RUN pip install --no-cache-dir --upgrade -r /app/requirements.txt

# Copy the entire server directory contents into the container
COPY server/ /app/

# Hugging Face Spaces strictly requires port 7860
EXPOSE 7860

# Run the FastAPI app using main.py inside the container root
# Assumes your FastAPI instance inside main.py is named 'app'
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "7860"]