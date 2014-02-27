# Re-create the climate data database.

rm -f climate.db
sqlite3 climate.db < make-db.sql

for country in $(cat countries.txt)
do
    for var in $(cat vars.txt)
    do
        grep -v year $country-$var.csv | python insert.py $country $var | sqlite3 climate.db
    done
done
