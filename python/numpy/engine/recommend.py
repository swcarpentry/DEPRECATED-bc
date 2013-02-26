'''
Original example from Toby Segaran: "Programming Collective Intelligence"
Altered by Richard T. Guy (2010)
'''

from math import sqrt
import numpy

EPS = 1.0e-9 # Never use == for floats.

raw_scores = {

  'Bhargan Basepair' : {
    'Jackson 1999' : 2.5,
    'Chen 2002' : 3.5,
    'Rollins and Khersau 2002' : 3.0,
    'El Awy 2005' : 3.5,
    'Chen 2008' : 2.5,
    'Falkirk et al 2006' : 3.0
  },

  'Fan Fullerene' : {
    'Jackson 1999' : 3.0,
    'Chen 2002' : 3.5,
    'Rollins and Khersau 2002' : 1.5,
    'El Awy 2005' : 5.0,
    'Falkirk et al 2006' : 3.0,
    'Chen 2008' : 3.5
  },

  'Helen Helmet' : {
    'Jackson 1999' : 2.5,
    'Chen 2002' : 3.0,
    'El Awy 2005' : 3.5,
    'Falkirk et al 2006' : 4.0
  },

  'Mehrdad Mapping' : {
    'Chen 2002' : 3.5,
    'Rollins and Khersau 2002' : 3.0,
    'Falkirk et al 2006' : 4.5,
    'El Awy 2005' : 4.0,
    'Chen 2008' : 2.5
  },

  'Miguel Monopole' : {
    'Jackson 1999' : 3.0,
    'Chen 2002' : 4.0,
    'Rollins and Khersau 2002' : 2.0,
    'El Awy 2005' : 3.0,
    'Falkirk et al 2006' : 3.0,
    'Chen 2008' : 2.0
  },

  'Gail Graphics' : {
    'Jackson 1999' : 3.0,
    'Chen 2002' : 4.0,
    'Falkirk et al 2006' : 3.0,
    'El Awy 2005' : 5.0,
    'Chen 2008' : 3.5
  },

  'Stephen Scanner' : {
    'Chen 2002' :4.5,
    'Chen 2008' :1.0,
    'El Awy 2005' :4.0
  }
}

def prep_data(all_scores):
  '''
  Turn {person : {title : score, ...} ...} into NumPy array.
  Each row is a person, each column is a paper title.
  Note that input data is sparse (does not contain all person X paper pairs).
  '''

  # Names of all people in alphabetical order.
  people = all_scores.keys()
  people.sort()

  # Names of all papers in alphabetical order.
  papers = set()
  for person in people:
    for title in all_scores[person].keys():
      papers.add(title)
  papers = list(papers)
  papers.sort()

  # Create and fill array.
  ratings = numpy.zeros((len(people), len(papers)))
  for (person_id, person) in enumerate(people):
    for (title_id, title) in enumerate(papers):
      rating = all_scores[person].get(title, 0)
      ratings[person_id, title_id] = float(rating)

  return people, papers, ratings

def sim_distance(prefs, left_index, right_index):
  '''
  Calculate distance-based similarity score for two people.
  Prefs is array[person X paper].

  Calculated a similarity difference btween two people (rows),
  which is 0 if they have no preferences in common.
  '''

  # Where do both people have preferences?
  left_has_prefs = prefs[left_index, :] > 0
  right_has_prefs = prefs[right_index, :] > 0
  mask = numpy.logical_and(left_has_prefs, right_has_prefs)

  # Not enough signal.
  if numpy.sum(mask) < EPS:
    return 0

  # Return sum-of-squares distance.
  diff = prefs[left_index, mask] - prefs[right_index, mask]
  sum_of_squares = numpy.linalg.norm(diff) ** 2
  result = 1. / (1. + sum_of_squares)
  return result

def sim_pearson(prefs, left_index, right_index):
  '''
  Calculate Pearson correlation between two individuals.
  '''

  # Where do both have ratings?
  rating_left = prefs[left_index, :]
  rating_right = prefs[right_index, :]
  mask = numpy.logical_and(rating_left > 0, rating_right > 0)

  # Note that summing over Booleans gives number of Trues
  num_common = sum(mask)

  # Return zero if there are no common ratings.
  if num_common == 0:
    return 0

  # Calculate Pearson score "r"
  varcovar = numpy.cov(rating_left[mask], rating_right[mask])
  numerator = varcovar[0,1]

  denominator = sqrt(varcovar[0,0]) * sqrt(varcovar[1,1])

  if denominator < EPS:
    return 0

  r = numerator / denominator
  return r

def top_matches(ratings, person, num, sim_func):
  '''
  Return the most similar individuals to a person.
  '''

  scores = []
  for other in range(ratings.shape[0]):
    if other != person:
      scores.append((sim_func(ratings, person, other), other))

  scores.sort()
  scores.reverse()
  return scores[0:num]

def calculate_similar(paper_ids, ratings, num=10):
  '''
  Find the papers that are most similar to each other.
  '''

  result = {}
  ratings_by_paper = ratings.T
  for item in range(ratings_by_paper.shape[0]):
    unnamed_scores = top_matches(ratings_by_paper, item, num, sim_distance)
    scores = [(x[0], paper_ids[x[1]]) for x in unnamed_scores]
    result[paper_ids[item]] = scores

  return result

def recommend(prefs, subject, sim_func):
  '''
  Get recommendations for an individual from a weighted average of other people.
  '''

  totals = {}
  sim_sums = {}
  num_people = prefs.shape[0]
  num_papers = prefs.shape[1]

  for other in range(num_people):

    # Don't compare people to themselves.
    if other == subject:
      continue
    sim = sim_func(prefs, subject, other)

    # ignore scores of zero or lower
    if sim < EPS:
      continue

    for title in range(num_papers):
      
      # Only score papers this person hasn't seen yet.
      if prefs[subject, title] < EPS and prefs[other, title] > EPS:
        
        # Similarity * Score
        if title in totals:
          totals[title] += prefs[other, title] * sim
        else:
          totals[title] = 0

        # Sum of similarities
        if title in sim_sums():
          sim_sums[title] += sim
        else:
          sim_sums[title] = 0

  # Create the normalized list
  
  rankings = []
  for title, total in totals.items():
    rankings.append((total/sim_sums[title], title))

  # Return the sorted list
  rankings.sort()
  rankings.reverse()
  return rankings

def test():
  person_ids, paper_ids, all_ratings = prep_data(raw_scores)
  print 'person_ids', person_ids
  print 'paper_ids', paper_ids
  print 'all_ratings', all_ratings
  print 'similarity distance', sim_distance(all_ratings, 0, 1)
  print 'similarity Pearson', sim_pearson(all_ratings, 0, 1)
  print top_matches(all_ratings, 0, 5, sim_pearson)
  print calculate_similar(paper_ids, all_ratings)
  print recommend(all_ratings, 0, sim_distance)
  print recommend(all_ratings, 1, sim_distance)
	
if __name__ == '__main__':
  test()
