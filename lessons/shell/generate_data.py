import sys
import os
import math
from random import choice, randint, random
import calendar

class Person:
    maxCI = 25
    # teenagers are hereby declared to be between 11 and 20 years old
    birthyears = range(1991,2000)
    repeatFraction = 0.1
    
    names = ['john', 'paul', 'george', 'ringo',\
        'baby','scary','posh','ginger','madonna',\
        'prince','robyn','beyonce','jay'] 
    words =['Beatle','Spice','Backstreet','Sync','Jonas',\
        'Lennon','McCartney','Starr','Harrison','Z',\
        'Carrot','Broccoli','Asparagus','Beet']
    CIs=range(1,maxCI+1)
    birthmonths= range(1,13)
    #ensure unique ids
    serialNum=173
    sexes=['M','F','N']

    def age(self, curyr=2011, curmo=11):
        return curyr+(1.*curmo-1.)/12. - self.birthyear - 1.*(self.birthmonth-1.)/12.

    def __init__(self):
        self.subject = choice(Person.names)+choice(Person.words)+ ('%03d' % Person.serialNum)
        Person.serialNum = Person.serialNum + 1

        self.birthyear  = choice(Person.birthyears)
        self.birthmonth = choice(Person.birthmonths)
   
        self.sex = choice(Person.sexes)
        age = self.age(2011,11)
        self.CI = choice(Person.CIs) 

        # newer CIs have better volume, discrimination;
        # range goes down with age.  (say). 

        CInewness = (self.CI-1.)/(1.*max(Person.CIs))
        # from oldest CI to newest, gain 2 volume pts: 
        self.trueVolume = randint(0,4)+randint(1,4)+round(2.*CInewness)

        # from oldest CI to newest, gain 3 discrimination pts: 
        self.trueDiscrimination = randint(0,3)+randint(1,4)+round(3.*CInewness)
        
        # 21-year-olds would lose 3 range points over 10 year olds (say)
        self.trueRange = randint(0,4)+randint(1,6)+round((10.-(self.age()-11.))*3./10.)

        # Most people don't repeat; those that do take the test 2-5 times
        if (random() > Person.repeatFraction):
            self.repeats = 1
        else:
            self.repeats=choice(range(2,6))


from numpy import polyfit, array
def test_peopleCorrelations():
    testpeople = []
    npeople = 4000
    for pnum in xrange(1,npeople):
        testpeople.append(Person())

    data = [[p.age(), p.CI, p.trueVolume, p.trueRange, p.trueDiscrimination] for p in testpeople]
    ages, cis, vols, ranges, discs = zip(*data)

    CIVolParam, dummy   = polyfit(cis, vols, 1) 
    CIRangeParam, dummy = polyfit(cis, ranges, 1) 
    CIDiscParam, dummy  = polyfit(cis, discs, 1)

    AgeVolParam, dummy   = polyfit(ages, vols, 1) 
    AgeRangeParam, dummy = polyfit(ages, ranges, 1) 
    AgeDiscParam, dummy  = polyfit(ages, discs, 1) 

    assert CIVolParam > 0.75*(2./25.) and CIVolParam < 1.25*(2./25.)
    assert CIDiscParam > 0.75*(3./25.) and CIDiscParam < 1.25*(3./25.)
    assert AgeRangeParam < 0.75*(-3./10.) and AgeRangeParam > 1.25*(-3./10.)

    zeroTol = 0.03
    assert abs(CIRangeParam) < zeroTol
    assert abs(AgeVolParam)  < zeroTol
    assert abs(AgeDiscParam) < zeroTol



class Measurement:
    incompleteFraction = 0.05
    serialNum = 211
    def randomDate(self):
        hrs = range(8,17)
        mins = range(1,60)
        secs = range(1,60)
        months = range(5,10)

        month = choice(months)
        monthname = calendar.month_abbr[month]
        day = choice(range(1,calendar.monthrange(2011, month)[1]))
        dayname = calendar.day_abbr[calendar.weekday(2011, month, day)]
        hr = choice(hrs)
        min = choice(mins)
        sec = choice(secs)
        
        datestring = '%s %s %d %02d:%02d:%02d %s' % (dayname, monthname, day, hr, min, sec, '2011')
        return [datestring, month, day, hr, min, sec]

    def limit(self,n):
        if n < 1 :
            n = 1
        if n > 10 :
            n = 10
        return n 

    def __init__(self, p):
        """Generate a result"""
        self.person = p
        self.datestring, self.month, self.day, self.hr, self.min, self.sec = self.randomDate();

        self.serialNum = Measurement.serialNum
        Measurement.serialNum = Measurement.serialNum + 1

        # +/- 1 random measurement error
        self.volume = self.person.trueVolume + choice([-1,0,0,0,+1])
        self.range  = self.person.trueRange + choice([-1,0,0,0,+1])
        self.discrimination  = self.person.trueDiscrimination + choice([-1,0,0,0,+1])

        self.volume = self.limit(self.volume)
        self.range = self.limit(self.range)
        self.discrimination = self.limit(self.discrimination)

        # before this date, things were being recorded 0..9 rather than 1..10
        fixmonth = 8
        fixday = 18
        fixhr = 10

        fixdate = fixmonth*10000 + fixday*100 + fixhr 
        checkdate = self.month*10000 + self.day*100 + self.hr 
        if checkdate < fixdate:
            self.volume = self.volume - 1
            self.range = self.range - 1
            self.discrimination = self.discrimination - 1
    
        if (random() < Measurement.incompleteFraction):
            self.discrimination = None
        

    def __str__(self):
        text = '# ' + '\n'
        text += "%s: %s\n" % ( 'Reported', self.datestring )
        text += "%s: %s\n" % ( 'Subject',  self.person.subject )
        text += "%s: %4d/%02d\n" % ( 'Year/month of birth', self.person.birthyear,  self.person.birthmonth )
        text += "%s: %s\n" % ( 'Sex', self.person.sex )
        text += "%s: %d\n" % ( 'CI type', self.person.CI )
        text += "%s: %d\n" % ( 'Volume', self.volume )
        text += "%s: %d\n" % ( 'Range', self.range )
        if self.discrimination is None :
            text += "%s: \n" % ( 'Discrimination' )
        else:
            text += "%s: %d\n" % ( 'Discrimination', self.discrimination )
    
        return text

class Datataker:
    names = ['angela', 'JamesD', 'jamesm', 'Frank_Richard',\
        'lab183','THOMAS','alexander','Beth','Lawrence',\
        'Toni', 'gerdal', 'Bert', 'Ernie', 'olivia', 'Leandra',\
        'sonya_p', 'h_jackson'] 
    filenamestyles = ['data_%d','Data%04d','%d','%04d','audioresult-%05d']
    suffixstyles = ['.dat','.txt','','','.DATA']
    tookNotesFraction = 0.5
    notes = ['Took data on Thursday and Friday until 4pm;\nAll day saturday.\n',\
             'Contact Janice about new calibration for data in August.\n',\
             'Submission of hours last week shows only 7 hours because \none was spent cleaning the lab.\n',\
             'Had some trouble accessing data submission form on Saturday,\nso fewer submissions then.\n',\
             'Third subject had real problems with the discrimiation test, so omitted.\n',\
             'Discrimination test seems kind of flaky - had to skip in several cases\n',\
             'Fuse blew midway through this weeks data taking,\nfewer results than last week.\n']
    notefilenames = ['notes.txt','NOTES','ReadMe','misc.txt','About']

    def __init__(self):
        self.name = choice(Datataker.names)
        Datataker.names.remove(self.name)
        self.filenameprefix = choice(Datataker.filenamestyles)
        self.filenamesuffix = choice(Datataker.suffixstyles)
        self.measures = []
        self.tookNotes = False
        if (random() < Datataker.tookNotesFraction) :
            self.tookNotes = True 
            self.notes = choice(Datataker.notes)
            self.noteFilename = choice(Datataker.notefilenames)

    def addmeasurement(self,measurement):
        self.measures.append(measurement)

    def write(self):
        os.mkdir(self.name)
        os.chdir(self.name)

        if (self.tookNotes):
            fname = self.noteFilename
            file = open(fname, 'w')
            file.write(self.notes)
            file.close()

        for m in self.measures:
            fname = self.filenameprefix % m.serialNum + self.filenamesuffix
            file = open(fname, 'w')
            file.write(str(m))
            file.close()
        os.chdir('..')
            
 
def main():
    #test_peopleCorrelations()

    npeople = 300 # should generate ~ .9*300 + 3.5*.1*300 ~ 375 files
    nfiles = 351

    people = []
    for pnum in range(npeople):
        people.append(Person())

    measurements = []
    for p in people:
        for m in range(p.repeats):
            measurements.append(Measurement(p))

    nexperimenters = 7
    experimenters = []
    for i in range(nexperimenters):
        experimenters.append(Datataker())

    for fnum in xrange(min(len(measurements), nfiles)):
        ex = choice(experimenters)
        ex.addmeasurement(measurements[fnum]) 

    os.mkdir('data')
    os.chdir('data')
    for ex in experimenters:
        ex.write()
    os.chdir('..')

if __name__=='__main__':
    sys.exit(main())

