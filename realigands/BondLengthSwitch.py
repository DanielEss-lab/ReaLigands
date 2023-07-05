def switch(metal, bonded, distance=0):
    if 30 >= metal >= 21 and bonded == 8:
        return distance <= 1.7
    elif 48 >= metal >= 39 and bonded == 8:
        return distance <= 1.91
    elif 80 >= metal >= 57 and bonded == 8:
        return distance <= 1.96

    elif 27 >= metal >= 21 and bonded == 7:
        return distance <= 1.675
    elif 30 >= metal >= 28 and bonded == 7:
        return distance <= 1.72
    elif 45 >= metal >= 39 and bonded == 7:
        return distance <= 1.92
    elif 48 >= metal >= 46 and bonded == 7:
        return distance <= 1.78
    elif 76 >= metal >= 57 and bonded == 7:
        return distance <= 1.87
    elif 80 >= metal >= 77 and bonded == 7:
        return distance <= 1.78

    elif 27 >= metal >= 21 and bonded == 6:
        return distance <= 1.84
    elif 30 >= metal >= 28 and bonded == 6:
        return distance <= 1.76
    elif 45 >= metal >= 39 and bonded == 6:
        return distance <= 1.94
    elif 48 >= metal >= 46 and bonded == 6:
        return distance <= 1.87
    elif 76 >= metal >= 57 and bonded == 6:
        return distance <= 1.89
    elif 80 >= metal >= 77 and bonded == 6:
        return distance <= 1.95
