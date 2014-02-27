for country in AFG ARG AUS BRA CAN CHN FRA KAZ NGA PER RUS THA TUR USA ZMB
do
    for var in tas pr
    do
        echo $country $var
        curl -o $country-$var.csv http://climatedataapi.worldbank.org/climateweb/rest/v1/country/cru/$var/year/$country.CSV
    done
done
