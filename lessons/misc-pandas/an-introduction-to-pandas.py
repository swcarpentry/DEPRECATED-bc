
## Exercise 1:
# How would we get the second to last date (EDT) in the dataset?
# Combine head() and tail()
last_two_dates = data.EDT.tail(2)
second_to_last_date = last_two_dates.head(1)
print second_to_last_date
  

## Exercise 2:
# What is the range of temperatures in the dataset?

hottest_temp = data.max_temp.max()  # Highest of the highs  
coldest_temp = data.min_temp.min()  # Lowest of the lows  
print "Temperature range:", hottest_temp - coldest_temp, "degrees F"  

# Temperature range: 105 degrees F


## Exercise 3:
  
#  Print out the cloud cover for each day in May.
#  *Hint: you can make datetime objects with the `datetime(year, month, day)` function*
  
datetime(2012, 5, 1)  # May 1st of 2012
  
data[datetime(2012, 5, 1):datetime(2012, 5, 31)].cloud_cover

## Exercise 4:
  
#  Was there any November rain?
  
d = datetime(2012, 1, 1)
d.strftime("%B")
  
november_rain = False
for date_idx, row in data.iterrows():
    if date_idx.strftime("%B") == "November" and "Rain" in row["events"]:
        november_rain = True
  
if november_rain:
      print "There was rain in November"
else:
      print "There was *not* rain in November"

 
  ## Exercise 5:
  
  
#  We'll replace "T" with a very small number, and convert the rest of the strings to floats:
  
# Convert precipitation to floating point number
# "T" means "trace of precipitation"
def precipitation_to_float(precip_str):
    if precip_str == "T":
        return 1e-10  # Very small value
    return float(precip_str)
  
data.precipitation = data.precipitation.apply(precipitation_to_float)
data.precipitation.head()


## Exercise 6:
  
#  Was the mean temperature more variable on days with rain and snow than on days with just rain or just snow?
  
days_with_rain = data[data.rain == True]
days_with_snow = data[data.snow == True]

rain_std = days_with_rain.mean_temp.std()
snow_std = days_with_snow.mean_temp.std()
  
if rain_std > snow_std:
   print "Rainy days were more variable"
elif snow_std > rain_std:
   print "Snowy days were more variable"
else:
   print "They were the same"

  
## Exercise 7:
  
#  Add the mean temperature to the previous plot using a green line. Also, add a legend with the `legend()` method of `ax`.
  
ax = data.max_temp.plot(title="Min and Max Temperatures")
data.min_temp.plot(style="red", ax=ax)
data.mean_temp.plot(style="green", ax=ax)
ax.set_ylabel("Temperature (F)")
  


