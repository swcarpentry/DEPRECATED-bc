# Setting credential for Git (Windows OS)

> I do not want to get prompt for my username and password every single time I push my commit to Github! 

The following gives the instruction to set up credential for Git for Windows users using [Windows Credential Store for Git](http://gitcredentialstore.codeplex.com/). 

The following simply paraphrases the instruction in [Windows Credential Store for Git](http://gitcredentialstore.codeplex.com/).

## Instruction

1. Download the [git-credential-winstore.exe](http://gitcredentialstore.codeplex.com/releases/view/103679) application.

2. Run it! It should work if GIT is in your "PATH" environment variable. If not, go to the directory where you download the application and run the following:
  
  ```
  git-credential-winstore -i "C:\Path\To\Git.exe"
  ```
For example, the path can be: "C:\Program Files (x86)\Git\bin\git.exe".

3. For the first time, you will be prompted to enter the username and password of your Github account, but otherwise, you are ready to go.

## Removing credential 

You may go to Git bash and type the following:
```
git config --unset --global credential.helper
```
