from data import *
from util import *


if __name__ == "__main__":

    while True:

        settings["STS"] = None

        while type(settings["STS"]) != int:
            try:
                settings["STS"] = int(input(settings["Question"]))
            except ValueError:
                print(settings["warning"])

        if settings["STS"] == 0:

            q1 = float(input("Введите коэффицент x1: "))
            q2 = float(input("Введите коэффицент x2: "))
            analyze_single_condition(q1, q2, **data)

        elif settings["STS"] == 1:
            graph(-1, 4, 1e-2, **data)
            plt.show()
            if input("Хотите ли перейти к расчету конкретной пары x1, x2? (Y/N)").lower() == "y":
                q1 = float(input("Введите коэффицент x1: "))
                q2 = float(input("Введите коэффицент x2: "))
                analyze_single_condition(q1, q2, **data)

        elif settings["STS"] == 2:
            planet(**data)

        elif settings["STS"] == 3:
            exit()