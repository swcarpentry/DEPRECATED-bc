##Software Carpentry Instructors Training

###Some ideas for the sessions, tasks and exercises
(work in progress to be improved)  


####**1. Making screencasts**  
There should be a **pre-defined set of lessons for making screencasts**. They should all come from the Software Carpentry "core". The trainees could only pick a topic from this set.  
Reasons:  
a) helps the trainees make up their mind;  
b) it is about seeing oneself teach but also they are supposed to teach SWC so the topic *should* be from SWC;  
c) the screencasts may be compared (if the participants want to do it - it would be a voluntary step).

After the screencasts are done and the participants received the feedback from the other two from the group, they could be regrouped based on the choice of the topic they decided to teach. That is those who made a screencast on "Loops in Python" would get together and so on. They don't have to watch their videos (in fact there is no time for that) but discuss their experience and the feedback they received. They could also **compare their strategies and approaches to teaching the topic**.

The set of topics should be send a few days before and the participants should be reminded that they will have to record themselves.   
Reasons:  
a) in Toronto some people didn't have enough time to prepare overnight for making the screencast;  
b) some participants said they would have brought better recording devices.

####**2. Observing Greg teaching a selected module from the SWC core**
The trainees are students. Greg teachers a part of a selected module from the SWC core (for example, something on command line, "head" and "tail"). The **teaching itself would be 30-45 mins**. 3-4 trainees could be helpers (or the instructors present at the training would be helpers, Aleksandra, Rob Dave etc.).   
The trainees get to observe and learn from The Most Experienced and The Best Instructor :-)  **All trainees have feedback forms which they fill in**.  
Once the teaching is finished, there is a group discussion on what was good, what was bad, which teaching approaches and tricks the trainees picked up from Greg.

####**3. Multiple choice questions with wrong answers**
The goal of this task is to develop way of analysing of the potential misconceptions which students may have. The instructors need to also become aware how they themselves may induce such misconceptions.  
How this task could be implemented:  
a) in groups of 3 the trainees come up with 3 questions and 4 answers to them; and explanations for the wrong answers ("why the student could choose one of the wrong answers");    
b) the groups exchange the questions with answers between them (but not reveal the explanation). The other group has to come up with the explanations for the wrong answers. This exercise will allow the group that came up with the wrong answers to see whether they can "simulate" student's mind and understand possible misconceptions. Also, trying to understand what could be the explanation behind the wrong that the other group created is a more realistic exercise - as if they actually tried to guess the reason behind the wrong answer. 

####**4. Strategies for preparing to teach material that you're not familiar with**
Strategies:  
a) come clean - tell the hosts, the administator and the other instructor which module(s) you would rather not teach because you are still not very confident with the material;  
b) negotiate with the other instructor splitting the material so that each one of you teaches the bit that you know/are comfortable with;  
c) train yourself up! Go through the *whole* material of a given module and learn it; write down all bits that you don't understand. Then go through it again, as if you were teaching - this will reveal even more bits of material that you don't quite get! Are you really sure what happens in Git when you create a branch?  

	TIP: It will take you 3-4 times longer than you estimate to prepare for teaching material that you are not familiar with or that you have never taught before. 
		
d) there is a big difference between "I (kind of) understand how that works and can do it" and "I can explain to someone else how that works and how to do it". As an instructor you should be able to do the latter. Pick up an exercise from the material that you are trying to learn and try to come up with a different example. How difficult is it for you to think of a different exercise? How does it help the students to deal with possible misconceptions?

####**5. Being a lead instructor - what that means and how to do it**    
"It sounds more scary than it really is." The goal of this session is to make the trainees more comfortable with the responsibilities and the work of the lead instructors; and to make them ready for this role.  
The structure of the session:  
a)   trainees read through the section on the website on "how to be a lead instructor" (this could be a homework to be done in the evening of Day1 or even sent out earlier).   
b)  They write down things that they:  
		- don't understand / don't know how to do;  
		- think they will struggle doing them.  
c) Then they work in groups of 3-4. Each person picks one thing from the above two lists (so two issues) and the group tries to find answers/solutions to them.  
d) Greg or Aleksandra talk about the "mentorship" - getting help (via email) from more experienced instructors who were the leads several times.


####**6. Creating a bootcamp repo from the template**  
This task is to simulate what normally a lead instructor has to do. I know that the repo will change, that the structure will change. But the purpose of this exercise is to allow the trainees to see how long it takes, identify what problems they may face (technical - needing help with Git; contents - needing to liaise with hosts to get info about venue, hours etc.; ) They will also (hopefully) see that it's not that scary. And even if we change the whole structure, they will have more confidence they can do it. Also, they should be encouraged to provide feedback (PRs!) on the instructions written for setting up a repo, find bugs  etc.  
[maybe the trainees could also be asked to set up an Etherpad and link it to the website they create?]
[This could be homework after Day1. Or homework after the Train the Trainers event]

####**7. How NOT to teach**  
The goal of this session is to make the trainees aware what kinds of their behaviour, wording, style of teaching may cause more damage than good.  
a) Telling students that they are doing "rubbish work" because they use Excel, Word; they don't modularise their code; don't test it etc. etc.    
[TO DO pull out that email sent to the Discuss list with student's feedback who felt overly criticized and discouraged for using Excel]  
b) Perpetually complaining about Windows/Mac and praising Linux/Unix - juxtaposing the former as "for amateurs/not real OS/broken OS" with the latter ("for professionals/for *real* geeks/for real,serious tasks")  
c) Criticizing GUI appications and describing command-line tools as "the only ones that do the job". Even worse, an instructor putiing himself/herself in a position of a "guru": "I'm a real programmer; don't know anything about that silly Windows stuff." By doing this you intimidate the students. Many of them have already been intimidated and came to SWC hoping to get help.    
d) Not managing expectations. Start with the hosts when you negotiate the curriculum. Make sure that they understand the main goals and purpose of SWC (we target beginners; it's not about the intrisics of object oriented programming). It's usually the leading instructors job but the SWC admin may be able to help you, if the hosts don't seem to be getting the message through. It is really crucial to make it clear to them because it's them who will be advertising/recommenting SWC workshop to the students. If they misadvertise it, you may end up with an audience who will be very disappointed they are taught how to modularise their code whilst they actually wanted to learn some MPI ("I thought it was supposed to be about scientific programming!"). Once you sorted the hosts, you should still make it clear to the attendees at the beginning of Day 1 what you and the other instructor will be doing over the course of the following two days and why you do it this way.  
e) Diving into discussing really complex and difficult questions asked by 1-2 people in the audience who clearly don't really need SWC. They already know it all and want to discuss whether Scala is better than Haskell. If you, as an instructor, engage in these questions because you feel they are far more interesting than teaching how to create functions in Python, you probably should not be an SWC instructor. If teaching beginners is *boring* for you and you'd rather have advanced discussions, instructing at SWC workshops is probably not be best way you could contribute to SWC. Maybe you'd be interested in debugging a number of srcipts which help to run our infrastructure? Or help with the admintool?  
f) Hide your impostor syndrome. If you suffer from the impostor syndrome, learn how to use it to your advantage. Remember that SWC's main goal is to help people who need to improve relatively basic software engineering and development skills (though their research/domain knowledge is probably very good). They will feel more at ease if you aren't afraid to admit that you don't think you're an expert on a given topic. But don't overdo it - you don't want to give an impression that SWC instructors know nothing about computers. We at least know how to switch them on.

+things that are on the SWC website under the header "what NOT to do"

####**8. How to help SWC grow as an instructor**   
At the beginning and at the end of each workshop encourage the attendees and helpers to get engaged in Software Carpentry. Have a slide ready showing different ways to engage and the path: attendee -> helper -> instructor. Have the information ready about the instructors training and point people to the right website (teaching.software-carpentry.org). Make sure that you challenge the misconception that an SWC instructor must be an experienced software developer with a degree in Computer Science. Give examples of some of the instructors who are just like the members of the audience. Maybe you yourself are an example like this?

####**9. Code of conduct** 

####**10. Dealing with audiences with a wide knowledge/skills range**  
Strategies:  
a) get the smart ones to be helpers. Occassionally, in collaboration with the hosts (of they know the audience very very well), you may be able to pair people up so that the beginners sit together with a more advanced person.  
b) (see 7.d);  
c) 

  

