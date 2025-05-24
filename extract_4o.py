import os
from dotenv import load_dotenv
from openai import OpenAI
import json
import csv
from tqdm import tqdm
from collections import defaultdict
import sys
import time

# Debugging?
DEBUG = False

# Get relative path to the script directory
script_dir = os.path.dirname(os.path.realpath(__file__))
dir_path = os.path.relpath(script_dir, os.getcwd())

# Load the OpenAI API key from a .env file
dotenv = os.path.join(dir_path, 'dot.env')
load_dotenv(dotenv)
client = OpenAI()

# Directory containing the example abstract/extract pairs
examples_dir = os.path.join(dir_path, 'examples')
if DEBUG:
    abstracts_dir = os.path.join(dir_path, 'test_abstracts')
    extracts_dir = os.path.join(dir_path, 'test_extracted')
    output_tsv = os.path.join(dir_path, 'test_extracted_data.tsv')
else:
    abstracts_dir = os.path.join(dir_path, 'abstracts')
    extracts_dir = os.path.join(dir_path, 'extracted')
    output_tsv = os.path.join(dir_path, 'extracted_data.tsv')

# Define measure categories
MEASURE_CATEGORIES = {
    "ratio": [
        "odds ratio", "or",
        "risk ratio", "rr",
        "hazard ratio", "hr",
        "prevalence ratio", "pr",
        "incidence ratio", "ir", "irr", 
        "incidence rate ratio",
        "incidence rate",
        "standardized incidence ratio", "sir",
        "relative risk",
        "rate ratio",
        "mortality ratio",
        "prevalence odds ratio"
    ],
    "association": [
        "regression coefficient", 
        "correlation coefficient",
        "linear model coefficient", 
        "logistic regression coefficient",
        "statistical test", 
        "interaction test",
        "regression analysis", 
        "correlation test",
        "covariance", "variance explained",
        "explained variation", "proportion explained",
        "percentage explained", "percentage",
        "factor loading", "factor score",
        "test for trend", "trend test",
        "trend analysis",
        "spearman correlation",
        "pearson correlation",
        "beta coefficient",
        "regression slope"
    ],
    "difference": [
        "mean difference", 
        "risk difference",
        "group comparison", 
        "statistical test",
        "anova", "analysis of variance",
        "t-test", "t test",
        "chi-square test", "chi square",
        "proportion difference",
        "rate difference"
    ]
}

# Define the function and how it should be structured within the API call
tools = [
    {
        "type": "function",
        "function": {
            "name": "extract_from_abstract",
            "description": "Extract structured data from abstracts. Must include PMID, PublicationDate, Author, Journal, and Measures. Each measure in Measures must include Metric, MeasureType, MeasureCategory, and (PValue) and/or (CIUpper AND CILower).",
            "parameters": {
                "type": "object",
                "properties": {
                    "PMID": {"type": "string", "description": "The PubMed ID from the abstract"},
                    "PublicationDate": {"type": "string", "description": "The publication date from the abstract"},
                    "Author": {"type": "string", "description": "The first author of the publication"},
                    "Journal": {"type": "string", "description": "The journal in which the publication appeared"},
                    "StudyDesign": {"type": "string", "description": "The study design of the publication; e.g., randomized controlled trial, cohort study, case-control study, etc."},
                    "Measures": {
                        "type": "array",
                        "description": "An array of measurements with point estimates and confidence intervals",
                        "items": {
                            "type": "object",
                            "properties": {
                                "PValue": {"type": "string", "description": "The p-value associated with the measurement; e.g.,'0.05' or '<0.01'"},
                                "PointEstimate": {"type": "number", "description": "The point estimate from the abstract (numeric value)"},
                                "CILower": {"type": "number", "description": "The lower bound of the confidence interval (numeric value)"},
                                "CIUpper": {"type": "number", "description": "The upper bound of the confidence interval (numeric value)"},
                                "ConfidenceLevel": {"type": "string", "description": "Confidence level of the interval, e.g., 95%, 99%"},
                                "MeasureType": {"type": "string", "description": "Type of statistical measure (Odds Ratio, Hazard Ratio, Relative Risk, Regression Coefficient, etc.)"},
                                "MeasureCategory": {"type": "string", "description": "Category of measure: ratio, association, or difference"},
                                "Metric": {"type": "string", "description": "The measurement with the point estimate and confidence interval"},
                            },
                            "required": ["Metric", "MeasureType", "MeasureCategory"]
                        }
                    }
                },
                "required": ["PMID", "PublicationDate", "Author", "Journal", "Measures"],
            },
        },
    }
]

# Define the function that will be called by the API
def extract_from_abstract(PMID, PublicationDate, author, journal, study_design, Measures):
    return {
        "PMID": PMID,
        "Date": PublicationDate,
        "Author": author,
        "Journal": journal,
        "StudyDesign": study_design,
        "Measures": Measures
    }

# API call to 4o-mini with the function calling
message_template=[
        {"role": "system", "content": f"""Extract statistical measures from abstracts, categorizing each according to these measure categories:

{json.dumps(MEASURE_CATEGORIES, indent=2)}

Format the output as JSON with this EXACT structure and field names:
{{
    "PMID": "string",
    "PublicationDate": "string",
    "Author": "string (lead author's last name)",  
    "Journal": "string",
    "StudyDesign": "string (optional)",
    "Measures": [
        {{
            "Metric": "string (what is being measured)",
            "MeasureType": "string (e.g., Odds Ratio, Mean difference)",
            "MeasureCategory": "string (ratio, association, or difference)",
            "PointEstimate": "number (optional)",
            "CILower": "number (optional)",
            "CIUpper": "number (optional)", 
            "ConfidenceLevel": "string (optional, e.g., 95%)",
            "PValue": "string (optional, e.g., 0.05 or <0.01)"
        }}
    ]
}}

IMPORTANT FIELD NAME RULES:
1. Use "Author" (not "LeadAuthor" or "Lead Author")
2. Field names are case-sensitive
3. No spaces in field names
4. Field names must match exactly as shown above

STATISTICAL SIGNIFICANCE RULES:
Each measure MUST include at least ONE of:
1. P-value ("PValue")
2. Confidence interval (both "CILower" and "CIUpper")

Examples of VALID measures (✅ has p-value OR confidence interval):
✅ {{
    "Metric": "Association with outcome",
    "MeasureType": "Odds Ratio",
    "MeasureCategory": "ratio",
    "PValue": "0.03",
    "PointEstimate": 1.5
}}

✅ {{
    "Metric": "Risk of disease",
    "MeasureType": "Relative Risk",
    "MeasureCategory": "ratio",
    "PointEstimate": 2.1,
    "CILower": 1.3,
    "CIUpper": 3.4
}}

Examples of INVALID measures (❌ missing both p-value AND confidence interval):
❌ {{
    "Metric": "Participation rate",
    "MeasureType": "Percentage",
    "MeasureCategory": "difference",
    "PointEstimate": 98,
    "PValue": null,
    "CILower": null,
    "CIUpper": null
}}

❌ {{
    "Metric": "Mean age of participants",
    "MeasureType": "Mean difference",
    "MeasureCategory": "difference",
    "PointEstimate": 45.2
}}

Examples of INVALID field names that will cause errors:
❌ "LeadAuthor": "Smith J"      # Wrong - use "Author" instead
❌ "Lead Author": "Smith J"     # Wrong - no spaces allowed
❌ "author": "Smith J"          # Wrong - case sensitive
❌ "AUTHORS": "Smith J"         # Wrong - must be "Author"
✅ "Author": "Smith"           # Correct format

Each measure must include:
1. Metric: What is being measured
2. MeasureType: Type of measure matching the categories above
3. MeasureCategory: ratio/association/difference
4. Statistical significance: P-value and/or confidence intervals (CILower + CIUpper)

Measures may be empty if no valid measures are found. If a measure is invalid (missing required fields or missing both p-value AND confidence intervals), exclude it from the output."""},
    ]

# Build the few-shot learning prompt
print("Building few-shot learning prompt...")
example_files = os.listdir(examples_dir)
for file in example_files:
    if file.endswith('.txt'):  # Identify abstract text files
        abstract_path = os.path.join(examples_dir, file)
        extract_path = os.path.join(examples_dir, file.replace('abstract', 'extract').replace('.txt', '.json'))

        # Read the abstract text
        with open(abstract_path, 'r') as file:
            abstract_text = file.read()
        
        # Load the corresponding extraction JSON
        with open(extract_path, 'r') as file:
            extraction_data = json.load(file)

        # Append the abstract and extraction data to the messages prompt
        message_template.append({"role": "user", "content": abstract_text})
        message_template.append({"role": "assistant", "content": json.dumps(extraction_data, indent=2)})

def scan_files():
    """Scan directories and count abstract/extract pairs"""
    abstracts = defaultdict(list)
    extracts = defaultdict(list)
    total_abstracts = 0
    existing_pairs = 0
    
    # Scan abstracts
    for root, _, files in os.walk(abstracts_dir):
        for file in files:
            if file.endswith('.txt'):
                rel_path = os.path.relpath(os.path.join(root, file), abstracts_dir)
                abstracts[os.path.dirname(rel_path)].append(file)
                total_abstracts += 1
    
    # Scan existing extracts
    for root, _, files in os.walk(extracts_dir):
        for file in files:
            if file.endswith('.json'):
                rel_path = os.path.relpath(os.path.join(root, file), extracts_dir)
                extracts[os.path.dirname(rel_path)].append(file)
                # Check if corresponding abstract exists
                abstract_file = file[:-5] + '.txt'
                if abstract_file in abstracts[os.path.dirname(rel_path)]:
                    existing_pairs += 1
    
    return total_abstracts, existing_pairs

# Define the function to extract structured data from an abstract using ChatGPT
def extract_abstract(abstract_text, message_template):
    messages = message_template.copy()
    messages.append({"role": "user", "content": abstract_text})
    
    try:
        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=messages,
            temperature=0,
            max_tokens=2048,
            tools=tools,
            tool_choice={"type": "function", "function": {"name": "extract_from_abstract"}}  # Force function call
        )

        response_message = response.choices[0].message
        
        # Check finish reason
        finish_reason = response.choices[0].finish_reason
        if finish_reason == "length":
            print("Error: The conversation was too long for the context window.")
            return []
        elif finish_reason == "content_filter":
            print("Error: The content was filtered due to policy violations.")
            return []
        
        # Get tool calls
        tool_calls = response_message.tool_calls
        if not tool_calls:
            print("Error: No tool calls found in response")
            return []

        # Parse the output from OpenAI API
        rows = []

        # Extract information for each tool call
        for tool_call in tool_calls:
            if tool_call.function.name != "extract_from_abstract":
                print(f"Warning: Unexpected function call to {tool_call.function.name}")
                continue

            try:
                args_js = json.loads(tool_call.function.arguments)
            except json.JSONDecodeError:
                print("Error: Invalid JSON in response")
                print("Response should be JSON:")
                print(tool_call.function.arguments)
                continue

            # Validate required fields
            required_fields = ['PMID', 'PublicationDate', 'Author', 'Journal', 'Measures']
            missing_fields = [field for field in required_fields if field not in args_js]
            if missing_fields:
                print(f"Error: Missing required fields: {missing_fields}")
                continue

            # Validate Measures is an array (can be empty)
            if not isinstance(args_js['Measures'], list):
                print("Error: Measures must be an array")
                continue

            # Extract data if validation passes
            pmid = args_js['PMID']
            date = args_js['PublicationDate']
            author = args_js['Author']
            journal = args_js['Journal']
            study_design = args_js.get('StudyDesign', None)
            measures = []

            # Only validate measures if there are any
            if args_js['Measures']:
                for measure in args_js['Measures']:
                    # Check required measure fields
                    measure_required = ['Metric', 'MeasureType', 'MeasureCategory']
                    missing_measure_fields = [field for field in measure_required if field not in measure]
                    if missing_measure_fields:
                        print(f"Warning: Measure missing required fields: {missing_measure_fields}")
                        continue

                    # Check statistical significance (must have p-value or confidence interval)
                    has_pvalue = 'PValue' in measure
                    has_ci = ('CILower' in measure and 'CIUpper' in measure)
                    if not (has_pvalue or has_ci):
                        print("Warning: Measure missing statistical significance (p-value or CI)")
                        continue

                    measures.append({
                        "Metric": measure['Metric'],
                        "MeasureType": measure['MeasureType'],
                        "MeasureCategory": measure.get('MeasureCategory', None),
                        "PointEstimate": measure.get('PointEstimate', None),
                        "CILower": measure.get('CILower', None),
                        "CIUpper": measure.get('CIUpper', None),
                        "ConfidenceLevel": measure.get('ConfidenceLevel', None),
                        "PValue": measure.get('PValue', None)
                    })

            # Add extraction regardless of measures (empty measures list is allowed)
            rows.append({
                "PMID": pmid,
                "Date": date,
                "Author": author,
                "Journal": journal,
                "StudyDesign": study_design,
                "Measures": measures
            })

    except Exception as e:
        print(f"Error during extraction: {str(e)}")
        return []

    return rows

print(f"\nExtracting from '{abstracts_dir}' to '{extracts_dir}'")
os.makedirs(extracts_dir, exist_ok=True)

# Scan files first
print("Scanning directories...")
total_abstracts, existing_pairs = scan_files()
remaining_abstracts = total_abstracts - existing_pairs
print(f"Found {total_abstracts} abstracts, {existing_pairs} already extracted")
print(f"Remaining to process: {remaining_abstracts}\n")

if remaining_abstracts == 0:
    print("All abstracts have been processed")
    sys.exit(0)

# Get list of remaining abstracts to process
remaining_files = []
for root, _, files in os.walk(abstracts_dir):
    for file in files:
        if not file.endswith('.txt'):
            continue
            
        abstract_path = os.path.join(root, file)
        extract_rel_path = os.path.relpath(abstract_path, abstracts_dir)
        extract_path = os.path.join(extracts_dir, extract_rel_path[:-4] + '.json')

        if not os.path.exists(extract_path):
            remaining_files.append((abstract_path, extract_path, extract_rel_path))

# Process remaining abstracts with tqdm progress bar
extraction_times = []
with tqdm(total=remaining_abstracts, desc="Extracting", unit="abstract") as pbar:
    for abstract_path, extract_path, extract_rel_path in remaining_files:
        os.makedirs(os.path.dirname(extract_path), exist_ok=True)

        extract_start = time.time()
        with open(abstract_path, 'r') as f:
            abstract_text = f.read()
        extracted_data = extract_abstract(abstract_text, message_template)

        with open(extract_path, 'w') as f:
            json.dump(extracted_data, f, indent=4)
        
        extract_time = time.time() - extract_start
        extraction_times.append(extract_time)
        
        # Update progress and show extraction details
        avg_time = sum(extraction_times) / len(extraction_times)
        pbar.set_postfix({"Last": f"{extract_time:.1f}s", "Avg": f"{avg_time:.1f}s"})
        pbar.update(1)

print("\nConverting to TSV...")
headers = ['PMID', 'Date', 'Author', 'Journal', 'Study_Design', 'Measure_Category', 'Point_Estimate', 'CI_Lower', 'CI_Upper', 'P_Value', 'Confidence_Level', 'Measure_Type', 'Metric']
with open(output_tsv, 'w', newline='', encoding='utf-8') as tsvfile:
    writer = csv.DictWriter(tsvfile, fieldnames=headers, delimiter='\t')
    writer.writeheader()
    
    # Walk through each file in the extracts directory
    for root, dirs, files in os.walk(extracts_dir):
        for filename in files:
            if filename.endswith('.json'):
                filepath = os.path.join(root, filename)
                with open(filepath, 'r') as f:
                    data = json.load(f)
                    for item in data:
                        for measure in item['Measures']:
                            writer.writerow({
                                'PMID': item['PMID'],
                                'Date': item['Date'],
                                'Author': item['Author'],
                                'Journal': item['Journal'],
                                'Study_Design': item.get('StudyDesign', None),
                                'Measure_Category': measure.get('MeasureCategory', None),
                                'Point_Estimate': measure.get('PointEstimate', None),
                                'CI_Lower': measure.get('CILower', None),
                                'CI_Upper': measure.get('CIUpper', None),
                                'P_Value': measure.get('PValue', None),
                                'Confidence_Level': measure.get('ConfidenceLevel', None),
                                'Measure_Type': measure['MeasureType'],
                                'Metric': measure['Metric']
                            })

print(f"Data conversion to TSV completed. Check {output_tsv} for the output.")
