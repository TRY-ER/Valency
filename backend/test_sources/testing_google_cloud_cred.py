from google import genai

API_KEY = None

# read .env to get the api key
with open("../.env") as f:
    lines = f.readlines()
    for l in lines:
        vals = l.strip().split("=")
        if len(vals) == 2:
            if vals[0] == "GOOGLE_CLOUD_PROJECT_API_KEY":
                API_KEY = vals[1]
                break

if API_KEY:
    client = genai.Client(api_key=API_KEY)
    response = client.models.generate_content(
        model="gemini-2.0-flash-exp", contents="Explain how AI works"
    )

    print("response >>", response.text)
else:
    print("could not get the api key")