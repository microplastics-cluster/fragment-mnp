import pandas as pd

data = pd.read_excel(r"Polymer Information for Model (1).xlsx")
df = pd.DataFrame(data)
d = {}

def dict_Initializer():
    for i in df["Polymer"]:
        if i not in d:
            d[i] = []
    for l in range(len(df["Polymer"])):
        check = list(df.loc[l, :])
        d[check[0]] = check[1:]

if __name__ == "__main__":
    dict_Initializer()
    
    for k,v in d.items():
        print("Key: %s and it's values %s\n" % (k, v))

