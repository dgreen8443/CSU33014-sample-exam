//
// CSU33014 Summer 2020 Additional Assignment
// Part B of a two-part assignment
//
// Please write your solution in this file

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include "csu33014-annual-partB-person.h"

void find_reachable_recursive(struct person *current, int steps_remaining,
                              bool *reachable)
{
  // mark current root person as reachable

  reachable[person_get_index(current)] = true;
  // now deal with this person's acquaintances
  if (steps_remaining > 0)
  {
    int num_known = person_get_num_known(current);
    for (int i = 0; i < num_known; i++)
    {
      struct person *acquaintance = person_get_acquaintance(current, i);
      find_reachable_recursive(acquaintance, steps_remaining - 1, reachable);
    }
  }
}

// computes the number of people within k degrees of the start person
int number_within_k_degrees(struct person *start, int total_people, int k)
{
  bool *reachable;
  int count;

  // maintain a boolean flag for each person indicating if they are visited
  reachable = malloc(sizeof(bool) * total_people);
  for (int i = 0; i < total_people; i++)
  {
    reachable[i] = false;
  }

  // now search for all people who are reachable with k steps
  find_reachable_recursive(start, k, reachable);

  // all visited people are marked reachable, so count them
  count = 0;
  for (int i = 0; i < total_people; i++)
  {
    if (reachable[i] == true)
    {
      count++;
    }
  }
  return count;
}
int BFS(struct person *pers, char *person_index, unsigned int k)
{
  struct list *checked = new_list();
  int level = 0;
  struct queue *queue = new_queue();
  add(pers, queue);
  add(NULL, queue);
  while (!isEmpty(queue) && level < k)
  {
    pers = poll(queue);
    if (pers == NULL)
    {
      printf("person null\n");
      level++;
      add(NULL, queue);
      if (peek(queue) == NULL)
        break;
      else
        continue;
    }
    printf("%s%s\n", "known..", pers->known_people[0]->person_index);
    for (int i = 0; i < pers->number_of_known_people; i++)
    {
      if (strcmp(pers->known_people[i]->person_index, person_index) == 0)
      {
        printf("%s%s%s\n", pers->person_index, " knows: ", person_index);
        return 1;
      }
      else
      {
        if (!is_in_list(pers->known_people[i]->person_index, checked))
        {
          add_to_list(pers->known_people[i]->person_index, checked);
          add(pers->known_people[i], queue);
        }
      }
    }
  }
  return 0;
}

int within_k_degrees(struct person *people[], int num_people, unsigned int k)
{
  for (int i = 0; i < num_people - 1; i++)
  {
    for (int j = i + 1; j < num_people; j++)
    {
      if (BFS(people[i], people[j]->person_index, k) == 0)
      {
        return 0;
      }
    }
  }
  return 1;
}

// computes the number of people within k degrees of the start person;
// less repeated computation than the simple original version
int less_redundant_number_within_k_degrees(struct person *start,
                                           int total_people, int k)
{
  //return number_within_k_degrees(start, total_people, k);
 within_k_degrees()

  

}

// computes the number of people within k degrees of the start person;
// parallel version of the code
int parallel_number_within_k_degrees(struct person *start,
                                     int total_people, int k)
{
  return number_within_k_degrees(start, total_people, k);
}
void find_reachable_recursive_parallel(struct person *current, int steps_remaining,
                                       bool *reachable)
{
  // mark current root person as reachable

  reachable[person_get_index(current)] = true;
  // now deal with this person's acquaintances
  if (steps_remaining > 0)
  {
    int num_known = person_get_num_known(current);
    for (int i = 0; i < num_known; i++)
    {
      struct person *acquaintance = person_get_acquaintance(current, i);
      find_reachable_recursive(acquaintance, steps_remaining - 1, reachable);
    }
  }
}
// computes the number of people within k degrees of the start person
int number_within_k_degrees_parallel(struct person *start, int total_people, int k)
{
  bool *reachable;
  int count;

  // maintain a boolean flag for each person indicating if they are visited
  reachable = malloc(sizeof(bool) * total_people);
  for (int i = 0; i < total_people; i++)
  {
    reachable[i] = false;
  }

  // now search for all people who are reachable with k steps
  find_reachable_recursive_parallel(start, k, reachable);

  // all visited people are marked reachable, so count them
  count = 0;
  for (int i = 0; i < total_people; i++)
  {
    if (reachable[i] == true)
    {
      count++;
    }
  }
  return count;
}
