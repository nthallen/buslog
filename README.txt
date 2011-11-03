
route tag title color oppositeColor latMin latMax lonMin lonMax
  stop tag title lat lon stopId (no idea what stopId is used for, but it's very similar to tag.)
  direction tag title name useForUI
    stop tag
  path
    point lat lon

Things to do:
  Update buslog to take a duration argument (hours) after which to quit
  Create shell script wrapper to create a subdirectory, run buslog, then
  rename and move the log file and delete the directory unless there is
  a "keep" file

Plot paths against streets to see what those are all about
Compare paths to direction
  Do I need to use path for proper distance calculations?
  Test whether distance to origination point is monotonically increasing
