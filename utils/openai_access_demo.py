import os
import time
import base64
from openai import OpenAI
import requests
import re
from utils.prompts import DEMO_SYSTEM_PROMPT

def get_oai_completion(prompt, stream=False, files=None):
    """
    Gets a completion from an OpenAI-compatible API.
    Handles both text and multimodal (image) inputs.
    """
    # Check if any files are uploaded
    if files and len(files) > 0:
        # Use the Qwen-VL model for image processing
        api_pools = [
            ("your_api_key","base_url","model_name")
        ]
        api = api_pools[0]
        api_key, base_url, model = api
        client = OpenAI(api_key=api_key, base_url=base_url)
        
        # Prepare the message content
        messages = [
            {"role": "system", "content": [{"type": "text", "text": DEMO_SYSTEM_PROMPT}]}
        ]
        
        user_content = []
        
        # Add images
        for file in files:
            print(file)
            if isinstance(file, str):  # It's a file path
                with open(file, "rb") as f:
                    image_data = f.read()
            else:  # Assume it's a file-like object
                image_data = file.read()
            
            # Convert to base64
            base64_image = base64.b64encode(image_data).decode('utf-8')
            user_content.append({
                "type": "image_url",
                "image_url": {
                    "url": f"data:image/jpeg;base64,{base64_image}"
                }
            })
        
        # Add the text prompt
        user_content.append({"type": "text", "text": prompt})
        
        messages.append({
            "role": "user",
            "content": user_content
        })
        
        try:
            response = client.chat.completions.create(
                model=model,
                messages=messages,
                temperature=0.8,
                max_tokens=8000,
                stream=stream
            )
            
            if stream:
                return response
            else:
                print(response.choices[0].message.content)
                return response.choices[0].message.content
                
        except Exception as e:
            print(f"Qwen-VL API error: {e}")
            return None
    else:
        # Original logic for text processing
        api_pools = [
            ("your_api_key","base_url","model_name")
        ]
        api = api_pools[0]
        api_key, base_url, model = api
        
        client = OpenAI(api_key=api_key, base_url=base_url)
        system_prompt = DEMO_SYSTEM_PROMPT
        user_prompt = prompt

        try:
            response = client.chat.completions.create(
                model=model,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt},
                ],
                temperature=0.8,
                max_tokens=8000,
                stream=stream
            )
            
            if stream:
                return response
            else:
                res = response.choices[0].message.content
                return res
        except Exception as e:
            print(e)
            return None

def call_chatgpt(ins, stream=False, files=None):
    """
    Calls the completion function with a retry mechanism.
    """
    success = False
    re_try_count = 5
    ans = ''
    while not success and re_try_count >= 0:
        re_try_count -= 1
        try:
            ans = get_oai_completion(ins, stream=stream, files=files)
            success = True
        except Exception as e:
            print(f"Retry times left: {re_try_count}; Error: {e}", flush=True)
            time.sleep(5)
    return ans

def call_chatgpt_stream(ins):
    """Calls the chat model in streaming mode."""
    return call_chatgpt(ins, stream=True)