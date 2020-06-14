from pymongo import MongoClient
import matplotlib.pyplot as plt
import pandas as pd

class DB_retriever:
    database_name = "porous"
    results_collection_name = "comp_results"

    def __init__(self):
        self._client = MongoClient()
        self._db = self._client[self.database_name]
        self._comp_res_col = self._db[self.results_collection_name]

    def extract_data(self, params):
        docs =  self._comp_res_col.find(params)
        Ras = []
        ampl_psis = []
        for doc in docs:
            Ras.append(doc["Ra"])
            ampl_psis.append(doc["ampl_psi"])
        return Ras, ampl_psis
    
    def extract_dataframe(self, params):
        docs =  self._comp_res_col.find(params, {"_id":0, "Ra":1, "Le":1, "m":1, "soret":1, "ampl":1, "omega":1, "ampl_psi":1})
        df = pd.DataFrame(docs)
        return df

if __name__ == "__main__":
    db = DB_retriever()
    params = {}
    print(db.extract_dataframe(params))