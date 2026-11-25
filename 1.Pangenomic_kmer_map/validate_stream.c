/* Validate and forward a TSV stream without a second disk pass. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

int main(int argc, char **argv) {
    if (argc != 4) return 2;
    uint64_t expected = strtoull(argv[1],NULL,10), cols = strtoull(argv[2],NULL,10);
    const char *header = argv[3]; size_t hlen = strlen(header), hpos = 0;
    static unsigned char buf[1024*1024];
    uint64_t rows=0, fields=1; int prev_delim=1; size_t n;
    while ((n=fread(buf,1,sizeof(buf),stdin))>0) {
        for (size_t i=0;i<n;i++) {
            unsigned char ch=buf[i];
            if (rows==0) {
                if (ch=='\n') { if(hpos!=hlen) goto bad; }
                else { if(hpos>=hlen || ch!=(unsigned char)header[hpos++]) goto bad; }
            }
            if (ch=='\t') { if(prev_delim) goto bad; fields++; prev_delim=1; }
            else if(ch=='\n') { if(prev_delim || fields!=cols) goto bad; rows++; fields=1; prev_delim=1; }
            else prev_delim=0;
        }
        if(fwrite(buf,1,n,stdout)!=n) return 3;
    }
    if(ferror(stdin) || rows!=expected || !prev_delim || fields!=1 || fflush(stdout)) goto bad;
    fprintf(stderr,"Validated %" PRIu64 " rows and %" PRIu64 " columns\n",rows,cols);
    return 0;
bad:
    fprintf(stderr,"Invalid TSV at or after row %" PRIu64 " (expected %" PRIu64 " rows)\n",rows+1,expected);
    return 1;
}
