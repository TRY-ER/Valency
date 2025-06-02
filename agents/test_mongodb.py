#!/usr/bin/env python3
"""
Test script for MongoDB integration with tool response storage and retrieval.
This script tests the MongoDBHandler class directly without going through the API.
"""

import os
import sys
import json
from dotenv import load_dotenv
from custom_utils.mongodb_handler import MongoDBHandler

# Load environment variables
load_dotenv()

def test_mongodb_connection():
    """Test MongoDB connection."""
    print("Testing MongoDB connection...")
    
    # Initialize MongoDB handler
    handler = MongoDBHandler()
    
    if not handler._initialized:
        print("❌ Failed to connect to MongoDB")
        return False
        
    print("✅ Successfully connected to MongoDB")
    return True

def test_store_and_retrieve():
    """Test storing and retrieving a tool response."""
    print("\nTesting store and retrieve functionality...")
    
    # Initialize MongoDB handler
    handler = MongoDBHandler()
    
    if not handler._initialized:
        print("❌ MongoDB connection not initialized, skipping test")
        return False
    
    # Generate a test tool ID
    test_tool_id = handler.generate_tool_id()
    print(f"Generated test tool ID: {test_tool_id}")
    
    # Create a test response
    test_response = {
        "isError": False,
        "content": [
            {"type": "text", "text": "This is a test response with very long content that would normally be truncated."}
        ]
    }
    
    # Store the test response
    mongo_id = handler.store_tool_response(
        tool_id=test_tool_id,
        tool_name="TestTool",
        user_id="test_user_123",
        session_id="test_session_456",
        response_data=test_response
    )
    
    if not mongo_id:
        print("❌ Failed to store test response")
        return False
        
    print(f"✅ Successfully stored test response with MongoDB ID: {mongo_id}")
    
    # Retrieve the test response
    retrieved = handler.get_tool_response(test_tool_id)
    
    if not retrieved:
        print("❌ Failed to retrieve test response")
        return False
    
    print("✅ Successfully retrieved test response:")
    print(f"Tool ID: {retrieved.get('tool_id')}")
    print(f"Tool Name: {retrieved.get('tool_name')}")
    print(f"User ID: {retrieved.get('user_id')}")
    print(f"Session ID: {retrieved.get('session_id')}")
    print(f"Created At: {retrieved.get('created_at')}")
    print(f"Response Data: {json.dumps(retrieved.get('response_data'), indent=2)}")
    
    return True

def test_get_session_responses():
    """Test retrieving all tool responses for a session."""
    print("\nTesting session responses retrieval...")
    
    # Initialize MongoDB handler
    handler = MongoDBHandler()
    
    if not handler._initialized:
        print("❌ MongoDB connection not initialized, skipping test")
        return False
    
    # Generate a unique session ID for this test
    test_session_id = f"test_session_{handler.generate_tool_id()}"
    print(f"Generated test session ID: {test_session_id}")
    
    # Store multiple test responses for the same session
    for i in range(3):
        test_tool_id = handler.generate_tool_id()
        test_response = {
            "isError": False,
            "content": [
                {"type": "text", "text": f"This is test response #{i+1} for session {test_session_id}"}
            ]
        }
        
        mongo_id = handler.store_tool_response(
            tool_id=test_tool_id,
            tool_name=f"TestTool{i+1}",
            user_id="test_user_123",
            session_id=test_session_id,
            response_data=test_response
        )
        
        if not mongo_id:
            print(f"❌ Failed to store test response #{i+1}")
            return False
            
        print(f"✅ Successfully stored test response #{i+1} with tool ID: {test_tool_id}")
    
    # Retrieve all responses for the session
    session_responses = handler.get_tool_responses_by_session(test_session_id)
    
    if not session_responses:
        print("❌ Failed to retrieve session responses")
        return False
    
    print(f"✅ Successfully retrieved {len(session_responses)} responses for session {test_session_id}")
    
    for i, response in enumerate(session_responses):
        print(f"\nResponse #{i+1}:")
        print(f"Tool ID: {response.get('tool_id')}")
        print(f"Tool Name: {response.get('tool_name')}")
        print(f"Created At: {response.get('created_at')}")
    
    return True

def main():
    """Run all tests."""
    print("=== MongoDB Integration Tests ===\n")
    
    # Test MongoDB connection
    if not test_mongodb_connection():
        print("\n❌ MongoDB connection test failed, exiting")
        sys.exit(1)
    
    # Test store and retrieve functionality
    if not test_store_and_retrieve():
        print("\n❌ Store and retrieve test failed")
    
    # Test session responses retrieval
    if not test_get_session_responses():
        print("\n❌ Session responses test failed")
    
    print("\n=== Tests completed ===")

if __name__ == "__main__":
    main()
