from pymongo import MongoClient

class DB_writer:
    database_name = "porous"
    results_collection_name = "comp_results"

    def __init__(self):
        self._client = MongoClient()
        self._db = self._client[self.database_name]
        self._comp_results_collection = self._db[self.results_collection_name]

    def insert_document(self, data):
        # TODO: check if document exists, drop document
        self._comp_results_collection.insert_one(data)
    
    def show_content(self):
        data = self._comp_results_collection.find()
        for p in data:
            print(p)
    
    def remove_collection(self):
        self._comp_results_collection.drop()

if __name__ == "__main__":
    db = DB_writer()
    #db.remove_collection()
    db.show_content()

