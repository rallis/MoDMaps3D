#cs
SEARCH
http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nuccore&id=[version]
RETRIEVE
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=[genid]&rettype=[type]&retmode=text
#ce

; BEGIN OF SETTINGS
$searchURL="http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nuccore&id="
$getURL="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="
$saveFormat=".fasta"
$fileName="human_chrom_dna"
$dirName=$fileName
; END OF SETTINGS


$file = FileOpen($fileName&".txt", 0)

If $file = -1 Then
    MsgBox(0, "Error", "Unable to open file.")
    Exit
EndIf

if Not FileExists($dirName) then
   DirCreate($dirName)
EndIf

$s=0
While 1
   $line = FileReadLine($file)
   If @error = -1 Then ExitLoop
	  $s=$s+1
   ;MsgBox(0,$s,$line)
   ConsoleWrite($s&": Downloading "&$line&@CRLF)
   InetGet($getURL&$line, $dirName&"\"&$line&$saveFormat,1,0)   
   
WEnd

FileClose($file)