-- Create a database to store climate values.

create table metadata(
    what    text    not null, -- name of measure
    units   text    not null, -- explanation of units
    note    text    not null  -- any other information
);

insert into metadata values('tas', 'deg. C', 'annual average surface temperature');
insert into metadata values('pr',  'mm',     'annual total precipitation');

create table readings(
    country text    not null, -- 3-letter ISO country code
    what    text    not null, -- name of observation
    year    integer not null, -- year of observation
    value   real    not null, -- recorded value
    foreign key(what) references metadata(what)
);
