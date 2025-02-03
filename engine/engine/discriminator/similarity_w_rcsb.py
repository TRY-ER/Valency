import requests
import json
from dataclasses import dataclass
from typing import List, Dict, Any
from urllib.parse import urlencode

RCSB_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

@dataclass
class RCSBQuery:
    entry_id: str
    assembly_id: str = "1"
    rows: int = 25
    start: int = 0

class RCSBSearcher:
    def __init__(self):
        self.base_url = RCSB_URL

    def build_query(self, query: RCSBQuery) -> Dict[str, Any]:
        return {
            "query": {
                "type": "terminal",
                "service": "structure",
                "parameters": {
                    "operator": "strict_shape_match",
                    "target_search_space": "assembly",
                    "value": {
                        "entry_id": query.entry_id,
                        "assembly_id": query.assembly_id
                    }
                }
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                    "start": query.start,
                    "rows": query.rows
                },
                "results_content_type": ["experimental"],
                "sort": [{"sort_by": "score", "direction": "desc"}],
                "scoring_strategy": "combined"
            }
        }

    def search(self, query: RCSBQuery) -> Dict[str, Any]:
        query_dict = self.build_query(query)
        params = {"json": json.dumps(query_dict)}
        url = f"{self.base_url}?{urlencode(params)}"
        
        response = requests.get(url)
        response.raise_for_status()
        return response.json()

# Usage example:
if __name__ == "__main__":
    searcher = RCSBSearcher()
    query = RCSBQuery(entry_id="1MO8", rows=1)
    results = searcher.search(query)
    print(json.dumps(results, indent=2))