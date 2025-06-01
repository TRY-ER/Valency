import httpx
from fastapi import HTTPException, Security
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
import os  # Import os to access environment variables
from dotenv import load_dotenv  # Import dotenv

load_dotenv()  # Load environment variables from .env file

# Use an environment variable for the auth service URL, with a default fallback
AUTH_SERVICE_URL = os.getenv("AUTH_SERVICE_URL", "http://localhost:8001/auth/users/me")
security = HTTPBearer()

async def verify_token(credentials: HTTPAuthorizationCredentials = Security(security)):
    token = credentials.credentials
    async with httpx.AsyncClient() as client:
        try:
            response = await client.get(AUTH_SERVICE_URL, headers={"Authorization": f"Bearer {token}"})
            response.raise_for_status()  # Raises an exception for 4XX/5XX responses
            return response.json()  # Or whatever user information you want to return
        except httpx.HTTPStatusError as e:
            raise HTTPException(status_code=e.response.status_code, detail="Could not validate credentials with auth service")
        except httpx.RequestError:
            raise HTTPException(status_code=503, detail="Auth service is unavailable")

