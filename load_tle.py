import os
import requests
from datetime import datetime

#ID_number must be an int which is the sat id number.
#25544 for ISS, for example.
#acceptable_age tells how old a tle can be before we choose to re-fetch
def get_tle(ID_number, acceptable_age = 3):
    ID_number = str(ID_number)
    if not os.path.exists(ID_number + ".tle"):
        return web_retrieve_tle(ID_number)
    else:
        with open(ID_number + ".tle") as f:
            loaded_tle = f.readlines()
            loaded_tle = [line[:-1] for line in loaded_tle]
        current_year = datetime.utcnow().timetuple().tm_year
        current_day_of_year = datetime.utcnow().timetuple().tm_yday
        current_epoch_day = str(current_year) + str(current_day_of_year)
        loaded_tle_epoch = "20" + loaded_tle[0][18:23]
        print(loaded_tle)
        tle_age = float(current_epoch_day) - float(loaded_tle_epoch)
        if tle_age > acceptable_age:
            print("Re-fetching old TLE")
            return web_retrieve_tle(ID_number)
        return loaded_tle
def web_retrieve_tle(ID_number):
    ID_number = str(ID_number)
    session = requests.session()
    url = "https://www.celestrak.com/NORAD/elements/gp.php?CATNR=" + ID_number
    page = session.get(url)
    sat_tle = page.text[:-2].split("\r\n")[1:]
    with open(ID_number + ".tle","w") as f:
        f.write("\n".join(sat_tle))
    return sat_tle

if __name__ == "__main__":
    print(get_tle(25544))