from fire import Fire

TARGET = "CTTGATCAT"
with open("Vibrio_cholerae.txt", "r") as f:
    GENOME = f.read()

result = []
for ix in range(len(GENOME)):
    if GENOME[ix:ix+9] == TARGET:
        result.append(ix)


if __name__ == "__main__":
    print(result)
