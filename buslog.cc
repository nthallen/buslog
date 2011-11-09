#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include "curllog/curl_obj.h"
#include "nortlib.h"

static int done = 0;

const char *find_last_time(xmlNodePtr tree) {
  const char *lastTime = NULL;
  for ( ; tree != NULL; tree = tree->next ) {
    if ( tree->name != NULL && strcmp((const char *)tree->name, "lastTime") == 0 ) {
      lastTime = (const char *)xmlGetProp(tree, (const xmlChar *)"time" );
      return lastTime;
    }
    lastTime = find_last_time( tree->children );
    if ( lastTime != NULL ) return lastTime;
  }
  return NULL;
}

static const char *get_prop( xmlNodePtr tree, const char *prop ) {
  const char *pval = (const char *)xmlGetProp(tree, (const xmlChar *)prop);
  if (pval == NULL) pval = "n/a";
  return pval;
}

void write_log_int( FILE *fp, const char *lastTime, xmlNodePtr tree ) {
  for ( ; tree != NULL; tree = tree->next ) {
    if ( tree->name != NULL && strcmp((const char *)tree->name, "vehicle") == 0) {
      const char *id = get_prop( tree, "id" );
      const char *dirTag = get_prop(tree, "dirTag");
      const char *lat = get_prop(tree, "lat");
      const char *lon = get_prop(tree, "lon");
      const char *secsSinceReport = get_prop(tree, "secsSinceReport");
      const char *heading = get_prop(tree, "heading");
      fprintf( fp, "\"%s\", \"%s\", %s, %s, %s, %s, %s\n",
        id, dirTag, lastTime, secsSinceReport, lat, lon, heading );
    }
    write_log_int( fp, lastTime, tree->children );
  }
}

void write_log( const char *logfile, const char *lastTime, xmlNodePtr tree ) {
  FILE *fp = fopen( logfile, "a" );
  if ( fp == NULL )
    nl_error(4, "Could not append to log file %s", logfile );
  write_log_int( fp, lastTime, tree );
  fclose(fp);
}

/**
 * Want to add command line options to select route, verbosity.
 */
int main( int argc, char **argv ) {
  const char *lastTime = "0";
  const char *route = "77";
  const char *logfile = "locations.log";
  double hours = 20.;
  time_t end_time, now;
  
  if ( argc > 1 ) route = argv[1];
  if ( argc > 2 ) hours = atof(argv[2]);
  printf("Logging data for MBTA route %s for %.1lf hours\n", route, hours);
  fflush(stdout);
  
  end_time = time(NULL) + (time_t)(hours * 3600);
  curl_obj co;
  co.set_log_level( CT_LOG_NOTHING );
  while (!done) {
    xmlNodePtr tree;
    char urlbuf[256];
    snprintf( urlbuf, 256,
      "http://webservices.nextbus.com/service/publicXMLFeed?command=vehicleLocations&a=mbta&r=%s&t=%s",
      route, lastTime );
    co.transaction_start("Get bus data");
    co.set_url(urlbuf);
    co.parse_request(1);
    co.perform("Get Locations");
    co.transaction_end();
    tree = co.get_parse_tree();
    lastTime = find_last_time(tree);
    if (lastTime != NULL) {
      write_log( logfile, lastTime, tree );
    } else lastTime = "0";
    co.write_transaction_log();
    now = time(NULL);
    if ( now >= end_time) break;
    sleep(30);
  }
  return 0;
}
