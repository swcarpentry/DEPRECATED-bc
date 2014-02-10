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
	quant   text,
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
