#
#  Tcl implementation of ErrLogger
#
#  Dave Brown, LBNL 4/7/2003
# Copyright Information:
#	Copyright (C) 2003 Lawrence Berkeley National Lab
#
#
#
#  Error-logging related variables
#


global ErrLogSeverityLevels
global ErrLogCurrentSeverity
global ErrLogNames

#
#  only define the error levels once
#
if { ! [info exists ErrLogCurrentSeverity] } {
    #
    #  Set the default error level as 'Routine'
    #
    set ErrLogCurrentSeverity 2
    #
    #  setup the names
    #
    lappend ErrLogSeverityLevels debugging trace routine warning error fatal
    lappend ErrLogNames Debugging Trace Routine Warning Error Fatal
}
#
#  Define the procs if not already done
#

if { [expr [llength [info procs "ErrMsg"]] == 0] } {

#
# test whether errors of the given severity would be reported, given the current
# configuration 
#
    proc ErrLogging { severity } {
	global ErrLogSeverityLevels
	global ErrLogCurrentSeverity
	#
	#  First, translate the severity name to a level
	#
	set iseverity [lsearch -exact $ErrLogSeverityLevels $severity]
	#
	#  Only produce output if the severity is equal to or greater than that requested
	#
	return "[expr [expr $iseverity == -1]  || [expr $iseverity >= $ErrLogCurrentSeverity]]"
    }
#
# standard messaging proc
#
    proc ErrMsg { severity message } {
	global ErrLogSeverityLevels
	global ErrLogCurrentSeverity
	global ErrLogNames

	set iseverity [lsearch -exact $ErrLogSeverityLevels $severity]
	if { [expr $iseverity == -1]  ||
	     [expr $iseverity >= $ErrLogCurrentSeverity] } {
	    #
	    if { [expr $iseverity == -1] } {
		puts "Invalid severity level $severity: printing message"
		set iseverity 4
	    }
	    #
	    #  Find the calling proc's name.  Skip it if it's 'ErrMsg'
	    #
	    set procname [lindex [info level [expr [info level] -1] ] 0]
	    if { $procname == "ErrMsg" } { set procname ""}
	    #
	    #  If the caller is just ErrMsg, set the caller to the script name
	    #
	    set fullscript [info script]
	    set script [string range $fullscript [expr [string last "/" $fullscript ] +1] end]
	    puts "${script}\:${procname}\:[lindex $ErrLogNames $iseverity]\: $message"
	    #
	    #  Abort if the level is fatal
	    #
	    if { $severity == "fatal" } { 
		puts "${script}\:${procname}\:Fatal\: tcl interpreter will now exit"
		exitTcl
	    }
	}
	return
    }
#
#  change (or report) the error logging severity
#
    proc ErrLoggingLevel { { severity ""} } {
	global ErrLogSeverityLevels
	global ErrLogCurrentSeverity
	#
	#  Make sure this level exists
	#
	if { $severity != "" } {
	    set iseverity [lsearch -exact $ErrLogSeverityLevels $severity]
	    
	    if { [expr $iseverity != -1] } {
		set ErrLogCurrentSeverity $iseverity
	    } else {
		ErrMsg error "Unknown severity level $severity"
	    }
	}
	return "[lindex $ErrLogSeverityLevels $ErrLogCurrentSeverity]"
    }
}
