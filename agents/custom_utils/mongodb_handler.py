import os
from typing import Optional, Dict, Any, List
from pymongo import MongoClient
from bson.objectid import ObjectId
import datetime
import uuid

class MongoDBHandler:
    """
    Handler for MongoDB operations to store and retrieve tool responses.
    """
    _instance = None
    
    def __new__(cls):
        """Singleton pattern to ensure only one MongoDB connection is created."""
        if cls._instance is None:
            cls._instance = super(MongoDBHandler, cls).__new__(cls)
            cls._instance._initialized = False
        return cls._instance
    
    def __init__(self):
        """Initialize MongoDB connection if not already initialized."""
        if self._initialized:
            return
            
        # Get MongoDB connection string from environment variable or use default
        mongo_uri = os.getenv("MONGODB_URI", "mongodb://localhost:27017/")
        db_name = os.getenv("MONGODB_DB_NAME", "valency_agents")
        
        try:
            self.client = MongoClient(mongo_uri)
            self.db = self.client[db_name]
            self.tool_responses = self.db["tool_responses"]
            print(f"MongoDB connection established to {db_name}")
            self._initialized = True
        except Exception as e:
            print(f"Error connecting to MongoDB: {str(e)}")
            self._initialized = False
    
    def store_tool_response(self, tool_id: str, tool_name: str, user_id: str, 
                           session_id: str, response_data: Any) -> str | None:
        """
        Store a complete tool response in MongoDB.
        
        Args:
            tool_id (str): Unique identifier for the tool call
            tool_name (str): Name of the tool that was called
            user_id (str): ID of the user making the request
            session_id (str): ID of the session
            response_data (Any): The complete tool response data
            
        Returns:
            str: MongoDB document ID of the stored response
        """
        if not self._initialized:
            print("MongoDB connection not initialized")
            return None
            
        # Create a document to store
        document = {
            "tool_id": tool_id,
            "tool_name": tool_name,
            "user_id": user_id,
            "session_id": session_id,
            "response_data": response_data,
            "created_at": datetime.datetime.utcnow(),
        }
        
        try:
            result = self.tool_responses.insert_one(document)
            return str(result.inserted_id)
        except Exception as e:
            print(f"Error storing tool response: {str(e)}")
            return None
    
    def get_tool_response(self, tool_id: str) -> Optional[Dict[str, Any]]:
        """
        Retrieve a tool response from MongoDB by tool_id.
        
        Args:
            tool_id (str): Unique identifier for the tool call
            
        Returns:
            Optional[Dict[str, Any]]: The complete tool response document, or None if not found
        """
        if not self._initialized:
            print("MongoDB connection not initialized")
            return None
            
        try:
            result = self.tool_responses.find_one({"tool_id": tool_id})
            return result
        except Exception as e:
            print(f"Error retrieving tool response: {str(e)}")
            return None
    
    def get_tool_responses_by_session(self, session_id: str) -> List[Dict[str, Any]]:
        """
        Retrieve all tool responses for a specific session.
        
        Args:
            session_id (str): ID of the session
            
        Returns:
            List[Dict[str, Any]]: List of tool response documents for the session
        """
        if not self._initialized:
            print("MongoDB connection not initialized")
            return []
            
        try:
            results = list(self.tool_responses.find({"session_id": session_id}))
            return results
        except Exception as e:
            print(f"Error retrieving tool responses for session: {str(e)}")
            return []
    
    def delete_tool_responses_by_session(self, session_id: str) -> int:
        """
        Delete all tool responses for a specific session.
        
        Args:
            session_id (str): ID of the session to delete responses for
            
        Returns:
            int: Number of documents deleted
        """
        if not self._initialized:
            print("MongoDB connection not initialized")
            return 0
            
        try:
            result = self.tool_responses.delete_many({"session_id": session_id})
            print(f"Deleted {result.deleted_count} tool responses for session {session_id}")
            return result.deleted_count
        except Exception as e:
            print(f"Error deleting tool responses for session {session_id}: {str(e)}")
            return 0
    
    def generate_tool_id(self) -> str:
        """
        Generate a unique tool ID.
        
        Returns:
            str: A unique ID for a tool call
        """
        return str(uuid.uuid4())
