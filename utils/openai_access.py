import openai
import time
from openai import OpenAI
# from zhipuai import ZhipuAI
import requests
import re
# from utils.parse_llm_output import try_parse_json_object

def get_oai_completion(prompt):
    api_pools = [
        ("your_api_key","base_url","model_name")
    ]
    api = api_pools[0]
    api_key, base_url, model = api
    
    if "GLM" in model:
        
        client = ZhipuAI(api_key=api_key)
        # from utils.prompts import GLM_JSON_RESPONSE_PREFIX, GLM_JSON_RESPONSE_SUFFIX, system_prompt
        # system_prompt = f"{GLM_JSON_RESPONSE_PREFIX}{system_prompt}"
        # user_prompt = f"{prompt}{GLM_JSON_RESPONSE_SUFFIX}"
        system_prompt = "You are a helpful assistant."
        user_prompt = prompt
    else:
        client = OpenAI(api_key=api_key, base_url=base_url)
        system_prompt = "You are a helpful assistant."
        user_prompt = prompt

    try:
        if "GLM" in model:
            response = client.chat.completions.create(
                model=model,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt},
                ],
                # response_format={ "type": "json_object" },
                temperature=0.1,
                top_p=0.7,
                stream=False
            )
        else:
            # print(user_prompt)
            response = client.chat.completions.create(
                model=model,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt},
                ],
                # response_format={ "type": "json_object" },
                #根据任务的不同来调整
                temperature=0.8,
                max_tokens=8000,
                stream=False
            )
        res = response.choices[0].message.content

        # if "GLM" in model:
        #     pattern = re.compile(r"```(?:json\s+)?(\{.*?\})\s*```", re.DOTALL)
        #     match = pattern.search(res)
        #     if match:
        #         gpt_output, _ = try_parse_json_object(match.group(1).strip())
        #     else:
        #         gpt_output = res
        # else:
        #     gpt_output = res
            
        pattern = re.compile(r"```(?:json\s+)?(\{.*?\})\s*```", re.DOTALL)
        match = pattern.search(res)
        # if match:
        #     gpt_output, _ = try_parse_json_object(match.group(1).strip())
        # else:
        #     gpt_output = res
        gpt_output = res

        return gpt_output
        
    except requests.exceptions.Timeout:
        print("The API request timed out. Please try again later.")
        return None
    except Exception as e:
        print(e)
        return None

def call_chatgpt(ins):
    success = False
    re_try_count = 5
    ans = ''
    while not success and re_try_count >= 0:
        re_try_count -= 1
        try:
            ans = get_oai_completion(ins)
            success = True
        except Exception as e:
            print(f"Retry times: {re_try_count}; Error: {e}", flush=True)
            time.sleep(5)
    return ans