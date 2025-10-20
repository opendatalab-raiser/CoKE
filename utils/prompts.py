ENZYME_PROMPT = """**You are a senior systems biologist.** Analyze the input information to predict ec number using structured reasoning. Crucially, implement a **self-correct mechanism** with these steps:

### Self-Correct Protocol
1. **Enzyme Verification**  
   - Discard ANY information contradicting the enzyme nature (catalytic activity).  
   - Example: If a GO term implies non-enzymatic function (e.g., structural role), reject it immediately.

2. **Conflict Resolution (Majority Rule)**  
   - Identify conflicts between:  
     - Motif vs. Motif  
     - GO term vs. GO term  
     - Motif vs. GO term  
   - **Resolution Principle**:  
     - If one element (A) conflicts with ≥2 logically consistent elements (B,C,D), discard A.  
     - Preserve high-confidence information supported by multiple sources.  
   - *Note*: Compatible functions (e.g., catalytic activity + cofactor binding) are NOT conflicts.

3. **Output Filtered Information**  
   - Explicitly list retained/discarded items with reasons before analysis.

### Final Output Requirement for EC Number

After completing the full biological analysis, you **must** conclude your entire response with a special section for automated parsing. This section must adhere to the following precise logic and format:

**Decision Logic:**

1.  **Default to a Single EC Number:** Your primary goal is to predict the **single, most likely EC number** for the protein's primary catalytic activity.
2.  **Handling Ambiguity:** If the evidence suggests a single function but points to several possible EC numbers (e.g., a family motif describes related but distinct activities), you must **commit to one choice**. Select the EC number that is most representative, most common, or best supported by the combined evidence. **Do not list multiple options out of uncertainty.**
3.  **Exception for Bifunctionality:** You may only predict multiple EC numbers if there is **explicit and strong evidence that a single protein is bifunctional**, meaning it contains distinct domains that perform two or more separate catalytic reactions. This requires clear support, such as a motif description explicitly stating "bifunctional" or the presence of multiple, distinct top-level catalytic GO terms (e.g., both a kinase and a cyclase activity).

**Formatting Rules:**

1.  The section must begin on a new line with the exact tag: `[EC_PREDICTION]`
2.  **Single Prediction (Standard Case):** Follow the tag with a single space and the predicted EC number.
    *   Example: `[EC_PREDICTION] 1.14.99.54`
3.  **Bifunctional Prediction (Exceptional Case):** List the EC numbers separated by a comma with no spaces.
    *   Example: `[EC_PREDICTION] 2.7.1.1,4.6.1.1`
4.  Do not add any other text, explanation, or punctuation on this line.
"""

FUNCTION_PROMPT = """**You are a senior systems biologist.** Analyze the input information to answer the given question.
"""

LLM_SCORE_PROMPT = """As an expert biologist, you are assigned to check one paragraph is aligned with facts or not. You will receive some facts, and
one paragraph. Score the paragraph between 0 to 100.
The score should be the format of {"score": score}
Here's the facts:
{{ground_truth}}
Here's the paragraph:
{{llm_answer}}
"""





