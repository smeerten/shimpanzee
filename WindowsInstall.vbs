Set oWS = WScript.CreateObject("WScript.Shell") 
set fso = CreateObject("Scripting.FileSystemObject")
strDesktop= oWS.SpecialFolders("Desktop") 
useDir = fso.GetAbsolutePathName(".")
sLinkFile = strDesktop + "\shimpanzee.lnk"  
Set oLink = oWS.CreateShortcut(sLinkFile) 
Target = useDir + "\WindowsRun.bat"
oLink.TargetPath = """"& Target &"""" 
oLink.IconLocation = useDir + "\logo.ico"
oLink.WorkingDirectory = useDir
oLink.Save 

' StartMenu
strStartMenu= oWS.SpecialFolders("Programs")
StartLocation = strStartMenu+"\shimpanzee.lnk"
 
Set oLink = oWS.CreateShortcut(StartLocation) 
oLink.TargetPath = """"& Target &"""" 
oLink.IconLocation = useDir + "\logo.ico"
oLink.WorkingDirectory = useDir
oLink.Save 
