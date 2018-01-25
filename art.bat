@echo off

REM
REM Batch file for starting Artemis on windows
REM

REM execute Artemis
start javaw -classpath .;lib\biojava.jar;lib\jobcontrol.jar;lib\jemAlign.jar uk.ac.sanger.artemis.components.ArtemisMain

