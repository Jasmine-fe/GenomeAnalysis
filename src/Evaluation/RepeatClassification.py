from prettytable import PrettyTable


class RepeatClassification:
    def __init__(self):
        pass

    def printFamilyNameTable(self, familyNameCounter):
        table = PrettyTable(["FamilyName", "Count"])
        for i in list(familyNameCounter):
            table.add_row([i[0], i[1]])
        print(table)
