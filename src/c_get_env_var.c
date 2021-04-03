#include <stdlib.h>
#include <string.h>

void get_environment_variable_(char c[1], int *max_char, char var)
{
  int i;
  for (i=0; i<*max_char; ++i)
    c[i] = ' ';
  strcpy(c, getenv(var));
  c[strlen(getenv(var))] = ' ';
}

// void get_environment_variable1_(c, max_char)
// char c[1];
// int *max_char;
// {
//   int i;
//   for (i=0; i<=*max_char-1; ++i) c[i]=' ';    
//   strcpy(c,getenv("DATA_DIRECTORY")); 
//   c[strlen(getenv("DATA_DIRECTORY"))]=' ';  
// }

// void get_environment_variable2_(c, max_char)
// char c[1];
// int  *max_char;
// {
//   int i;
//   for (i=0; i<=*max_char-1; ++i) c[i]=' ';    
//   strcpy(c,getenv("RESTART_DIRECTORY")); 
//   c[strlen(getenv("RESTART_DIRECTORY"))]=' ';  
// }
