-- The `Person` table is used to explain the most basic queries.
-- Note that `danforth` has no measurements.
create table Person(
	ident    text,
	personal text,
	family	 text
);

insert into Person values('dyer',     'William',   'Dyer');
insert into Person values('pb',       'Frank',     'Pabodie');
insert into Person values('lake',     'Anderson',  'Lake');
insert into Person values('roe',      'Valentina', 'Roerich');
insert into Person values('danforth', 'Frank',     'Danforth');

-- The `Site` table is equally simple.  Use it to explain the
-- difference between databases and spreadsheets: in a spreadsheet,
-- the lat/long of measurements would probably be duplicated.
create table Site(
	name text,
	lat  real,
	long real
);

insert into Site values('DR-1', -49.85, -128.57);
insert into Site values('DR-3', -47.15, -126.72);
insert into Site values('MSK-4', -48.87, -123.40);

-- `Visited` is an enhanced `join` table: it connects to the lat/long
-- of specific measurements, and also provides their dates.
-- Note that #752 is missing a date; we use this to talk about NULL.
create table Visited(
	ident integer,
	site  text,
	dated text
);

insert into Visited values(619, 'DR-1',  '1927-02-08');
insert into Visited values(622, 'DR-1',  '1927-02-10');
insert into Visited values(734, 'DR-3',  '1939-01-07');
insert into Visited values(735, 'DR-3',  '1930-01-12');
insert into Visited values(751, 'DR-3',  '1930-02-26');
insert into Visited values(752, 'DR-3',  null);
insert into Visited values(837, 'MSK-4', '1932-01-14');
insert into Visited values(844, 'DR-1',  '1932-03-22');

-- The `Survey` table is the actual readings.  Join it with `Site` to
-- get lat/long, and with `Visited` to get dates (except for #752).
-- Note that Roerich's salinity measurements are an order of magnitude
-- too large (use this to talk about data cleanup).  Note also that
-- there are two cases where we don't know who took the measurement,
-- and that in most cases we don't have an entry (null or not) for the
-- temperature.
create table Survey(
	taken   integer,
	person  text,
	quant   real,
	reading real
);

insert into Survey values(619, 'dyer', 'rad',    9.82);
insert into Survey values(619, 'dyer', 'sal',    0.13);
insert into Survey values(622, 'dyer', 'rad',    7.80);
insert into Survey values(622, 'dyer', 'sal',    0.09);
insert into Survey values(734, 'pb',   'rad',    8.41);
insert into Survey values(734, 'lake', 'sal',    0.05);
insert into Survey values(734, 'pb',   'temp', -21.50);
insert into Survey values(735, 'pb',   'rad',    7.22);
insert into Survey values(735, null,   'sal',    0.06);
insert into Survey values(735, null,   'temp', -26.00);
insert into Survey values(751, 'pb',   'rad',    4.35);
insert into Survey values(751, 'pb',   'temp', -18.50);
insert into Survey values(751, 'lake', 'sal',    0.10);
insert into Survey values(752, 'lake', 'rad',    2.19);
insert into Survey values(752, 'lake', 'sal',    0.09);
insert into Survey values(752, 'lake', 'temp', -16.00);
insert into Survey values(752, 'roe',  'sal',   41.60);
insert into Survey values(837, 'lake', 'rad',    1.46);
insert into Survey values(837, 'lake', 'sal',    0.21);
insert into Survey values(837, 'roe',  'sal',   22.50);
insert into Survey values(844, 'roe',  'rad',   11.25);

select '========================================';
select 'Selecting';

select '----------------------------------------';
select 'get scientist names';
select family, personal from Person;

select '----------------------------------------';
select 'commands are case insensitive';
SeLeCt famILY, PERSonal frOM PERson;

select '----------------------------------------';
select 'we control column order';
select personal, family from Person;

select '----------------------------------------';
select 'repeat columns';
select ident, ident, ident from Person;

select '----------------------------------------';
select 'use * for wildcard';
select * from Person;

select '========================================';
select 'Removing Duplicates';

select '----------------------------------------';
select 'show quantities in survey table';
select quant from Survey;

select '----------------------------------------';
select 'unique quantity names';
select distinct quant from Survey;

select '----------------------------------------';
select 'tuple uniqueness';
select distinct taken, quant from Survey;

select '========================================';
select 'Filtering';

select '----------------------------------------';
select 'when a particular site was visited';
select * from Visited where site='DR-1';

select '----------------------------------------';
select 'filtering columns after "where"';
select ident from Visited where site='DR-1';

select '----------------------------------------';
select 'when a particular site was visited after 1930';
select * from Visited where (site='DR-1') and (dated>='1930-00-00');

select '----------------------------------------';
select 'using "or" instead of "and"';
select * from Survey where person='lake' or person='roe';

select '----------------------------------------';
select 'using "in" instead of "or"';
select * from Survey where person in ('lake', 'roe');

select '----------------------------------------';
select 'with parentheses';
select * from Survey where quant='sal' and person='lake' or person='roe';

select '----------------------------------------';
select 'using distinct with "in"';
select distinct person, quant from Survey where person='lake' or person='roe';

select '========================================';
select 'Calculating New Values';

select '----------------------------------------';
select 'correct radiation readings';
select 1.05 * reading from Survey where quant='rad';

select '----------------------------------------';
select 'convert temperatures to Celsius';
select taken, round(5*(reading-32)/9, 2) from Survey where quant='temp';

select '----------------------------------------';
select 'formatting names';
select personal || ' ' || family from Person;

select '========================================';
select 'Ordering Results';

select '----------------------------------------';
select 'ascending is the default';
select reading from Survey where quant='rad' order by reading;

select '----------------------------------------';
select 'order descending';
select reading from Survey where quant='rad' order by reading desc;

select '----------------------------------------';
select 'ordering and sub-ordering';
select taken, person from Survey order by taken asc, person desc;

select '----------------------------------------';
select 'removing duplicates';
select distinct taken, person from Survey order by taken asc, person desc;

select '----------------------------------------';
select 'sorted by undisplayed column';
select reading from Survey where quant='rad' order by taken;

select '----------------------------------------';
select 'randomization';
select random(), ident from Person;

select '----------------------------------------';
select 'random sorting';
select ident from Person order by random();

select '========================================';
select 'Missing Data';

select '----------------------------------------';
select 'all Visited data';
select * from Visited;

select '----------------------------------------';
select 'visits before 1930';
select * from Visited where dated<'1930-00-00';

select '----------------------------------------';
select 'visits after 1930';
select * from Visited where dated>='1930-00-00';

select '----------------------------------------';
select 'visits with unknown dates (wrong)';
select * from Visited where dated=NULL;

select '----------------------------------------';
select '!= does not work either';
select * from Visited where dated!=NULL;

select '----------------------------------------';
select 'visits with unknown dates (right)';
select * from Visited where dated is NULL;

select '----------------------------------------';
select 'visits with known dates';
select * from Visited where dated is not NULL;

select '----------------------------------------';
select 'salinity not measured by Lake';
select * from Survey where quant='sal' and person!='lake';

select '----------------------------------------';
select 'salinity not known to be measured by Lake';
select * from Survey where quant='sal' and (person!='lake' or person is null);

select '========================================';
select 'Aggregation';

select '----------------------------------------';
select 'all dates';
select dated from Visited;

select '----------------------------------------';
select 'date range';
select min(dated) from Visited;
select max(dated) from Visited;
select min(dated), max(dated) from Visited;

select '----------------------------------------';
select 'averaging';
select avg(reading) from Survey where quant='sal';

select '----------------------------------------';
select 'counting';
select count(reading) from Survey where quant='sal';

select '----------------------------------------';
select 'sum';
select sum(reading) from Survey where quant='sal';

select '----------------------------------------';
select 'averaging sensible values';
select min(reading), max(reading) from Survey where quant='sal' and reading<=1.0;

select '----------------------------------------';
select 'unaggregated with aggregated takes arbitrary';
select person, count(*) from Survey where quant='sal' and reading<=1.0;

select '----------------------------------------';
select 'what happens when there is no input?'
select person, max(reading), sum(reading) from Survey where quant='missing';

select '----------------------------------------';
select 'what happens when an input is null?'
select min(dated) from Visited;

select '----------------------------------------';
select 'what happens when an input is null?'
select min(dated) from Visited where dated is not null;

select '----------------------------------------';
select 'grouping Visited by site only keeps arbitrary';
select * from Visited group by site;

select '----------------------------------------';
select 'get date ranges for sites';
select   site, min(dated), max(dated) from Visited
group by site;

select '========================================';
select 'Grouping';

select '----------------------------------------';
select 'radiation reading for a particular person';
select person, count(reading), round(avg(reading), 2)
from  Survey
where quant='rad'
and   person='dyer';

select '----------------------------------------';
select 'radiation reading for an arbitrary person';
select person, count(reading), round(avg(reading), 2)
from  Survey
where quant='rad';

select '----------------------------------------';
select 'radiation readings by person';
select   person, count(reading), round(avg(reading), 2)
from     Survey
where    quant='rad'
group by person;

select '----------------------------------------';
select 'average by scientist and quantity';
select   person, quant, count(reading), round(avg(reading), 2)
from     Survey
group by person, quant;

select '----------------------------------------';
select 'grouping and ordering';
select   person, quant, count(reading), round(avg(reading), 2)
from     Survey
group by person, quant
order by person, quant;

select '----------------------------------------';
select 'grouping and ordering without null';
select   person, quant, count(reading), round(avg(reading), 2)
from     Survey
where    person is not null
group by person, quant
order by person, quant;

select '========================================';
select 'Combining Data';

select '----------------------------------------';
select 'combine "Site" with "Visited"';
select * from Site join Visited;

select '----------------------------------------';
select 'filter where sites match';
select * from Site join Visited on Site.name=Visited.site;

select '----------------------------------------';
select 'get latitude, longitude, and date';
select Site.lat, Site.long, Visited.dated
from   Site join Visited
on     Site.name=Visited.site;

select '----------------------------------------';
select 'get latitude, longitude, date, quantity, and reading';
select Site.lat, Site.long, Visited.dated, Survey.quant, Survey.reading
from   Site join Visited join Survey
on     Site.name=Visited.site
and    Visited.ident=Survey.taken
and    Visited.dated is not null;

select '----------------------------------------';
select 'get all radiation readings since 1930';
select 'but notice that #752 is missing (NULL)...';
select Survey.reading
from   Survey join Visited
where  Survey.taken=Visited.ident
  and  Survey.quant='rad'
  and  Visited.dated>='1930-00-00';

select '----------------------------------------';
select 'get row IDs from Person table';
select rowid, * from Person;
