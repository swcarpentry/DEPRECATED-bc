---
layout: lesson
root: ../..
title: Numbers
---
There are other ways to interact with remote files other than git.
In fact, if we wish to download files from the shell we can use
Wget, cURL, and ftp.

#### Wget

Wget is a simple tool developed for the GNU Project that downloads files with the HTTP, HTTPS and FTP protocols. It is widely used by unix users and is available with most Linux distributions.

To download this lesson from the web via http with the URL http://software-carpentry.org/v5/novice/extras/10-file_transfer.html we can simply type:

~~~
$ wget http://software-carpentry.org/v5/novice/extras/10-file_transfer.html
~~~
{:class="in"}

To produce output similar to:
~~~
--2014-11-21 09:41:31--  http://software-carpentry.org/v5/novice/extras/10-file_transfer.html
Resolving software-carpentry.org (software-carpentry.org)... 174.136.14.108
Connecting to software-carpentry.org (software-carpentry.org)|174.136.14.108|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 8901 (8.7K) [text/html]
Saving to: `10-file_transfer.html'

100%[======================================>] 8,901       --.-K/s   in 0.05s   

2014-11-21 09:41:31 (187 KB/s) - `10-file_transfer.html' saved [8901/8901]
~~~
{:class="out"}

### cURL

Alternatively, we can use cURL. It supports a much larger range of protocolsincluding common mail based protocols like pop3 and smtp. 

Functionally equivalent behaviour to the wget http request above may be invoked by typing:

~~~
$ curl -o 10-file_transfer.html http://software-carpentry.org/v5/novice/extras/10-file_transfer.html
~~~
{:class="in"}
To produce output similar to:
~~~
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                   Dload  Upload   Total   Spent    Left  Speed
                                   100 14005  100 14005    0     0  35170      0 --:--:-- --:--:-- --:--:--  105k
~~~
{:class="out"}

### lftp

Another option is lftp. It has a lot of capability, and even does simple bitorrent. 

~~~~
$ lftp -c get 03-review.html http://software-carpentry.org/v5/novice/extras/03-review.html
~~~
(:class="in"}
~~~
get: Not connected
~~~
(:class="out"}

Please refer to the man page to experiment with other protocols, such as FTP (e.g. man wget, man curl, man lftp).

Sources:

Wget - Wikipedia, http://en.wikipedia.org/wiki/Wget, Accessed: 11/21/2011

cURL - Wikipedia, http://en.wikipedia.org/wiki/CURL, Accessed: 11/21/2011

lftp - Wikipedia, http://en.wikipedia.org/wiki/Lftp, Accessed: 11/21/2011

lftp man page, Accessed: 11/21/2011



