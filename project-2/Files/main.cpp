#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <cstring>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

// constants
#define TRAINERS 200
#define USERS 100
#define MOVIES 1000
#define KSIMILAR 30
#define MAX_RATING 5

using namespace std;

// global vars
vector<vector<int> > training(TRAINERS, vector<int>(MOVIES, -1));
vector<vector<int> > users(USERS, vector<int>(MOVIES, -1));
vector<vector<int> > results(TRAINERS, vector<int>(MOVIES, -1));
int offset;
int testNum;

class Movie {
private:
    int _id;
    float _similarity;
    vector<int> _userRatings;
    int _reqRating;
    int _average;
public:
    Movie(int id, float similarity, vector<int> userRatings, int reqRating) {
        _id = id;
        _userRatings = userRatings;
        _similarity = similarity;
        _reqRating = reqRating;
        _average = calcAverage(userRatings);
    }

    // calculate average rating for movie
    int calcAverage(vector<int> userRatings) {
        int sum = 0;
        int count = 0;
        
        for(int i = 0; i < userRatings.size(); i++) {
            if(userRatings[i] > 0) {
                sum += userRatings[i];
                count++;
            }
            
        }

        if(sum <= 0)
            count = 1;
        
        return round(sum / count);
    }
    // getters
    float getSimilarity() { return _similarity; }
    int getReqRating() { return _reqRating; }
    float getAverage() { return _average; }
    int getRating(int i) { return _userRatings[i]; }

    // setters
    float setSimilarity(float sim) { _similarity = sim; }
};

class User {
private:
    int _reqRating;
    float _similarity;
    vector<int> _movieRatings;
    int _average;
public:
    // constructor
    User(vector<int> movieRatings, float similarity, int reqRating, vector<int> allMovieRatings) {
        _movieRatings = movieRatings;
        _similarity = similarity;
        _reqRating = reqRating;
        _average = calcAverage(allMovieRatings);
    }    

    // calculate average rating
    int calcAverage(vector<int> allMovieRatings) {
        int sum = 0;
        int count = 0;
        
        for(int i = 0; i < allMovieRatings.size(); i++) {
            if(allMovieRatings[i] > 0) {
                sum += allMovieRatings[i];
                count++;
            }
            
        }

        if(sum <= 0)
            count = 1;
        
        return round(sum / count);
    }

    // getters
    float getSimilarity() { return _similarity; }
    int getReqRating() { return _reqRating; }
    int getAverage() { return _average; }

    // setters
    float setSimilarity(float sim) { _similarity = sim; }
};

int itemBased(int user, int movie) {
    // user information
    vector<int> userMovieRatings;       // movie ratings from user for given movies 
    vector<int> userMovies;             // movies rated by user
    int rating = 0;                     // predicted rating

    // sim information
    vector<Movie> ratedMovies;
    vector<int> userRatings;
    vector<int> reqRatings;

    // similarity variables
    float dot = 0.0;
    float ru1i = 0.0;
    float ru2i = 0.0;
    float numerator = 0.0;
    float denominator = 0.0;

    // prediction variables
    float wau = 0.0;                    // similarity between movies
    int rui = 0;                        // similar movie's ratings

    // iterators
    int i, j, k, l;
    i = 0;
    j = 0;
    k = 0;
    l = 0;

    // find movies rated by user
    for(i = 0; i < MOVIES; i++) {
        if(users[user][i] != -1 && users[user][i] != 0) {                   // if user rated this movie
            userRatings.clear();                                            // clear userRatings vector
            userMovies.push_back(i);                                        // record movie #
            userMovieRatings.push_back(users[user][i]);                     // record user rating
            for(j = 0; j < TRAINERS; j++) {                                 // record all trainer's ratings for movie i
                userRatings.push_back(training[j][i]);
            }
            ratedMovies.push_back(Movie(i, 0.0, userRatings, users[user][i]));
        }
    }

    // create current User for average
    User curUser(userMovieRatings, 0.0, 0, userMovieRatings);

    // get all trainer ratings for current movie
    for(i = 0; i < TRAINERS; i++) {
        reqRatings.push_back(training[i][movie]);
    }

    // create current movie
    Movie curMovie(movie, 0.0, reqRatings, 0);

    // find similarity
    for(i = 0; i < testNum; i++) {

        wau = 0.0;
        dot = 0.0;
        ru1i = 0.0;
        ru2i = 0.0;

        for(j = 0; j < TRAINERS; j++) {                                                     // for all trainers
            if(curMovie.getRating(j) != 0 && ratedMovies[i].getRating(j) != 0) {            // if both movies were rated by the trainer
                dot += (float)((curMovie.getRating(j) - curUser.getAverage()) * (ratedMovies[i].getRating(j) - curUser.getAverage()));              // numer dot product
                ru1i += (float)pow((curMovie.getRating(j) - curUser.getAverage()), 2);            // denom left
                ru2i += (float)pow((ratedMovies[i].getRating(j) - curUser.getAverage()), 2);              // denom right
            }
        } 

        ru1i = sqrt(ru1i);
        ru2i = sqrt(ru2i);

        if(ru1i == 0 || ru2i == 0) {                    // if denom = 0
            wau = 0.0;                                  // then similarity = 0
        } else {
            wau = dot / (ru1i * ru2i);                  // else, compute similarity
        }

        ratedMovies[i].setSimilarity(wau);              // set similarity b/w current movie and similar movie

    }
    
    numerator = 0.0;
    denominator = 0.0;

    // user rating prediction calculation
    for(i = 0; i < testNum; i++) {                      // for x movies rated by the user
        wau = ratedMovies[i].getSimilarity();           // get movie's similarity
        rui = ratedMovies[i].getReqRating();            // get movie's needed rating
        numerator += (wau * (float)rui);                // calc numer
        denominator += abs(wau);                        // calc denom
    }

    if(denominator != 0)
        rating = round((numerator / denominator));

    if(rating <= 0) {                               // if no one else rated the movie
        rating = round(curUser.getAverage());       // rating = movie average
    }

    // restrict rating to a maximum
    if(rating > MAX_RATING) {
        rating = MAX_RATING;
    } 

    ratedMovies.clear();

    return rating;
}

int pearsonCase(int user, int movie) {
    // user information
    vector<int> userMovieRatings;       // movie ratings from user for given movies 
    vector<int> userMovies;             // movies rated by user
    int rating = 0;                     // predicted rating

    // sim information
    vector<int> simMovieRatings;                    // sim movie ratings for same movies as the user
    vector<int> allSimMovieRatings(MOVIES, 0);      // all sim movie ratings

    // relevant users (30)
    vector<User> relUsers;              // vector of relevant trainers
    int relUsersCount = 0;              // number of relevant trainers
    int reqRating = 0;                  // rating needed from trainer

    // similarity variables
    float dot = 0.0;
    float ru1i = 0.0;
    float ru2i = 0.0;
    float numerator = 0.0;
    float denominator = 0.0;

    // prediction variables
    float wau = 0.0;                    // similarity between user and trainer
    float wauPrime = 0.0;               // adjusted similarity between user and trainer
    float rho = 2.5;                    // case amplification power
    int rui = 0;                        // similar user's rating on movie

    // iterators
    int i, j, k, l;
    i = 0;
    j = 0;
    k = 0;
    l = 0;

    // find movies rated by user
    for(i = 0; i < MOVIES; i++) {
        if(users[user][i] != -1 && users[user][i] != 0) {               // if user rated this movie
            userMovies.push_back(i);                                    // record movie #
            userMovieRatings.push_back(users[user][i]);                 // record user rating
        }
    }

    // create current User
    User curUser(userMovieRatings, 0.0, 0, userMovieRatings);

    // find similar trainers
    for(i = 0; i < TRAINERS; i++) {                                     // for all trainers
        simMovieRatings.clear();                                        // clear movie ratings (same vector from previous trainer)
        for(j = 0; j < testNum; j++) {                                  // for all movies that were rated by the user
            if(training[i][userMovies[j]] < 1) {                        // if trainer didn't rate the movie
                simMovieRatings.push_back(0);                           // then set movie rating to 0
            } else {
                simMovieRatings.push_back(training[i][userMovies[j]]);        // get movie rating from trainer
            }
        }

        // get trainer's ratings
        for(j = 0; j < MOVIES; j++) {
            allSimMovieRatings[j] = training[i][j];
        }

        reqRating = training[i][movie];             // get trainer's rating 
        if(reqRating < 1) {                         // if trainer's rating is < 1, set to 0
            reqRating = 0;
        }

        // create similar user
        User simUser(simMovieRatings, 0.0, reqRating, allSimMovieRatings);

        dot = 0.0;
        ru1i = 0.0;
        ru2i = 0.0;

        // find similarity
        for(j = 0; j < testNum; j++) {
            if(userMovieRatings[j] != 0 && simMovieRatings[j] != 0) {               // if movie is rated by both users
                dot += (userMovieRatings[j] - curUser.getAverage()) * (simMovieRatings[j] - simUser.getAverage());                     // numer dot product
                ru1i += pow(userMovieRatings[j] - curUser.getAverage(), 2);                   // denom left
                ru2i += pow(simMovieRatings[j] - simUser.getAverage(), 2);                   // denom right
            }
        }
        ru1i = sqrt(ru1i);
        ru2i = sqrt(ru2i);

        if(ru1i == 0 || ru2i == 0) {                            // if denom = 0
            wau = 0.0;                                          // then similarity = 0
        } else {
            wau = dot / (ru1i * ru2i);                          // else, compute similarity
        }

        // case amplification
        wauPrime = wau * pow(abs(wau), rho);

        simUser.setSimilarity(wauPrime);             // set similarity b/w user and trainer

        if(reqRating != 0) {
            if(relUsersCount < KSIMILAR) {                  // if < KSIMILAR relevant trainers
                relUsers.push_back(simUser);                // record trainers ID
                relUsersCount++;                            // increment # of relevant trainers
            } else {                                        // else (already KSIMILAR trainers)
                for(k = 0; k < KSIMILAR; k++) {             // find least similar trainer
                    if(relUsers[k].getSimilarity() < relUsers[l].getSimilarity()) {
                        l = k;
                    }
                }
                relUsers[l] = simUser;                      // replace least similar trainer with current trainer
            }
        }
        
    }

    numerator = 0.0;
    denominator = 0.0;

    // user rating prediction calculation
    for(i = 0; i < relUsers.size(); i++) {                 // for k similar trainers
        wau = relUsers[i].getSimilarity();                     // get trainer's similarity
        rui = relUsers[i].getReqRating();                  // get trainer's needed rating
        numerator += wau * ((float)rui - relUsers[i].getAverage());              // calc numer
        denominator += abs(wau);                    // calc denom
    }

    if(denominator != 0)
        rating = curUser.getAverage() + round(numerator / denominator);

    if(rating <= 0) {                       // if no one else rated the movie
        rating = curUser.getAverage();      // rating = user average
    }

    // restrict rating to a maximum
    if(rating > MAX_RATING) {
        rating = MAX_RATING;
    } 

    relUsers.clear();

    return rating;
}

int pearsonInverse(int user, int movie) {
    // user information
    vector<int> userMovieRatings;       // movie ratings from user for given movies 
    vector<int> userMovies;             // movies rated by user
    int rating = 0;                     // predicted rating

    // sim information
    vector<int> simMovieRatings;                    // sim movie ratings for same movies as the user
    vector<int> allSimMovieRatings(MOVIES, 0);      // all sim movie ratings

    // relevant users (30)
    vector<User> relUsers;              // vector of relevant trainers
    int relUsersCount = 0;              // number of relevant trainers
    int reqRating = 0;                  // rating needed from trainer

    // similarity variables
    float iuf = 0.0;
    int m = TRAINERS;
    int mj = 1;
    float dot = 0.0;
    float ru1i = 0.0;
    float ru2i = 0.0;
    float numerator = 0.0;
    float denominator = 0.0;

    // prediction variables
    float wau = 0.0;                    // similarity between user and trainer
    int rui = 0;                        // similar user's rating on movie

    // iterators
    int i, j, k, l;
    i = 0;
    j = 0;
    k = 0;
    l = 0;

    // find movies rated by user
    for(i = 0; i < MOVIES; i++) {
        if(users[user][i] != -1 && users[user][i] != 0) {               // if user rated this movie
            userMovies.push_back(i);                                        // record movie #
            userMovieRatings.push_back(users[user][i]);                     // record user rating
        }
    }

    // create current User
    User curUser(userMovieRatings, 0.0, 0, userMovieRatings);

    // find similar trainers
    for(i = 0; i < TRAINERS; i++) {                     // for all trainers
        simMovieRatings.clear();
        for(j = 0; j < testNum; j++) {                      // for all movies that were rated by the user
            if(training[i][userMovies[j]] < 1) {                // if trainer didn't rate the movie
                simMovieRatings.push_back(0);                             // then set movie rating to 0
            } else {
                simMovieRatings.push_back(training[i][userMovies[j]]);        // get movie rating from trainer
            }
        }

        // get trainer's ratings
        for(j = 0; j < MOVIES; j++) {
            allSimMovieRatings[j] = training[i][j];
        }

        reqRating = training[i][movie];                 // get trainer's rating 
        if(reqRating < 1) {                 // if trainer's rating is < 1 (0 or -1)
            reqRating = 0;
        }

        // create similar user
        User simUser(simMovieRatings, 0.0, reqRating, allSimMovieRatings);

        dot = 0.0;
        ru1i = 0.0;
        ru2i = 0.0;

        // find similarity
        for(j = 0; j < testNum; j++) {
            if(userMovieRatings[k] != 0 && simMovieRatings[k] != 0) {               // if movie is rated by both users
                // find mj, number of users m who rated movie j
                mj = 0;
                for(k = 0; k < TRAINERS; k++) {
                    if(training[k][userMovies[j]] > 0) {
                        mj++;
                    }
                }
                
                // find iuf
                if(mj != 0) {
                    iuf = log2((float)m / (float)mj);
                } else {
                    iuf = 0.0;
                }

                dot += (((float)userMovieRatings[j]) - (float)curUser.getAverage()) * (((float)simMovieRatings[j] * iuf) - (float)simUser.getAverage());                 // numer dot product
                ru1i += pow(((float)simMovieRatings[j] * iuf) - (float)simUser.getAverage(), 2);                // denom left
                ru2i += pow(((float)userMovieRatings[j]) - (float)curUser.getAverage(), 2);               // denom right
            }
        }

        ru1i = sqrt(ru1i);
        ru2i = sqrt(ru2i);

        if(ru1i == 0 || ru2i == 0) {                            // if denom = 0
            wau = 0.0;                                          // then similarity = 0
        } else {
            wau = dot / (ru1i * ru2i);                          // else, compute similarity
        }

        simUser.setSimilarity(wau);             // set similarity b/w user and trainer

        if(reqRating != 0) {
            if(relUsersCount < KSIMILAR) {                  // if < KSIMILAR relevant trainers
                relUsers.push_back(simUser);                // record trainers ID
                relUsersCount++;                            // increment # of relevant trainers
            } else {                                        // else (already KSIMILAR trainers)
                for(k = 0; k < KSIMILAR; k++) {             // find least similar trainer
                    if(relUsers[k].getSimilarity() < relUsers[l].getSimilarity()) {
                        l = k;
                    }
                }
                relUsers[l] = simUser;                      // replace least similar trainer with current trainer
            }
        }
        
    }

    numerator = 0.0;
    denominator = 0.0;

    // user rating prediction calculation
    for(i = 0; i < relUsers.size(); i++) {                 // for k similar trainers
        wau = relUsers[i].getSimilarity();                     // get trainer's similarity
        rui = relUsers[i].getReqRating();                  // get trainer's needed rating
        numerator += wau * ((float)rui - relUsers[i].getAverage());              // calc numer
        denominator += abs(wau);                    // calc denom
    }

    if(denominator != 0)
        rating = curUser.getAverage() + round((numerator / denominator));

    if(rating <= 0) {                               // if no one else rated the movie
        rating = curUser.getAverage();       // rating = user average
    }

    // restrict rating to a maximum
    if(rating > MAX_RATING) {
        rating = MAX_RATING;
    } 

    relUsers.clear();

    return rating;
}

int pearsonCorrelation(int user, int movie) {
    // user information
    vector<int> userMovieRatings;       // movie ratings from user for given movies 
    vector<int> userMovies;             // movies rated by user
    int rating = 0;                     // predicted rating

    // sim information
    vector<int> simMovieRatings;                    // sim movie ratings for same movies as the user
    vector<int> allSimMovieRatings(MOVIES, 0);      // all sim movie ratings

    // relevant users (30)
    vector<User> relUsers;              // vector of relevant trainers
    int relUsersCount = 0;              // number of relevant trainers
    int reqRating = 0;                  // rating needed from trainer

    // similarity variables
    float dot = 0.0;
    float ru1i = 0.0;
    float ru2i = 0.0;
    float numerator = 0.0;
    float denominator = 0.0;

    // prediction variables
    float wau = 0.0;                    // similarity between user and trainer
    int rui = 0;                        // similar user's rating on movie

    // iterators
    int i, j, k, l;
    i = 0;
    j = 0;
    k = 0;
    l = 0;

    // find movies rated by user
    for(i = 0; i < MOVIES; i++) {
        if(users[user][i] != -1 && users[user][i] != 0) {               // if user rated this movie
            userMovies.push_back(i);                                    // record movie #
            userMovieRatings.push_back(users[user][i]);                 // record user rating
        }
    }

    // create current User
    User curUser(userMovieRatings, 0.0, 0, userMovieRatings);

    // find similar trainers
    for(i = 0; i < TRAINERS; i++) {                                     // for all trainers
        simMovieRatings.clear();                                        // clear movie ratings (same vector from previous trainer)
        for(j = 0; j < testNum; j++) {                                  // for all movies that were rated by the user
            if(training[i][userMovies[j]] < 1) {                        // if trainer didn't rate the movie
                simMovieRatings.push_back(0);                           // then set movie rating to 0
            } else {
                simMovieRatings.push_back(training[i][userMovies[j]]);        // get movie rating from trainer
            }
        }

        // get trainer's ratings
        for(j = 0; j < MOVIES; j++) {
            allSimMovieRatings[j] = training[i][j];
        }

        reqRating = training[i][movie];             // get trainer's rating 
        if(reqRating < 1) {                         // if trainer's rating is < 1, set to 0
            reqRating = 0;
        }

        // create similar user
        User simUser(simMovieRatings, 0.0, reqRating, allSimMovieRatings);

        dot = 0.0;
        ru1i = 0.0;
        ru2i = 0.0;

        // find similarity
        for(j = 0; j < testNum; j++) {
            if(userMovieRatings[j] != 0 && simMovieRatings[j] != 0) {               // if movie is rated by both users
                dot += (userMovieRatings[j] - curUser.getAverage()) * (simMovieRatings[j] - simUser.getAverage());                     // numer dot product
                ru1i += pow(userMovieRatings[j] - curUser.getAverage(), 2);                   // denom left
                ru2i += pow(simMovieRatings[j] - simUser.getAverage(), 2);                   // denom right
            }
        }
        ru1i = sqrt(ru1i);
        ru2i = sqrt(ru2i);

        if(ru1i == 0 || ru2i == 0) {                            // if denom = 0
            wau = 0.0;                                          // then similarity = 0
        } else {
            wau = dot / (ru1i * ru2i);                          // else, compute similarity
        }

        simUser.setSimilarity(wau);             // set similarity b/w user and trainer

        if(reqRating != 0) {
            if(relUsersCount < KSIMILAR) {                  // if < KSIMILAR relevant trainers
                relUsers.push_back(simUser);                // record trainers ID
                relUsersCount++;                            // increment # of relevant trainers
            } else {                                        // else (already KSIMILAR trainers)
                for(k = 0; k < KSIMILAR; k++) {             // find least similar trainer
                    if(relUsers[k].getSimilarity() < relUsers[l].getSimilarity()) {
                        l = k;
                    }
                }
                relUsers[l] = simUser;                      // replace least similar trainer with current trainer
            }
        }
        
    }

    numerator = 0.0;
    denominator = 0.0;

    // user rating prediction calculation
    for(i = 0; i < relUsers.size(); i++) {                 // for k similar trainers
        wau = relUsers[i].getSimilarity();                     // get trainer's similarity
        rui = relUsers[i].getReqRating();                  // get trainer's needed rating
        numerator += wau * (float)(rui - relUsers[i].getAverage());              // calc numer
        denominator += abs(wau);                    // calc denom
    }

    if(denominator != 0)
        rating = curUser.getAverage() + round((numerator / denominator));

    if(rating <= 0) {                       // if no one else rated the movie
        rating = curUser.getAverage();      // rating = user average
    }

    // restrict rating to a maximum
    if(rating > MAX_RATING) {
        rating = MAX_RATING;
    } 

    relUsers.clear();

    return rating;
}

int cosineSimilarity (int user, int movie)
{
    // user information
    vector<int> userMovieRatings;       // movie ratings from user for given movies 
    vector<int> userMovies;             // movies rated by user
    int rating = 0;                     // predicted rating

    // sim information
    vector<int> simMovieRatings;                    // sim movie ratings for same movies as the user
    vector<int> allSimMovieRatings(MOVIES, 0);      // all sim movie ratings

    // relevant users (30)
    vector<User> relUsers;              // vector of relevant trainers
    int relUsersCount = 0;              // number of relevant trainers
    int reqRating = 0;                  // rating needed from trainer

    // similarity variables
    float dot = 0.0;
    float ru1i = 0.0;
    float ru2i = 0.0;
    float numerator = 0.0;
    float denominator = 0.0;

    // prediction variables
    float wau = 0.0;                    // similarity between user and trainer
    int rui = 0;                        // similar user's rating on movie

    // iterators
    int i, j, k, l;
    i = 0;
    j = 0;
    k = 0;
    l = 0;

    // find movies rated by user
    for(i = 0; i < MOVIES; i++) {
        if(users[user][i] != -1 && users[user][i] != 0) {                   // if user rated this movie
            userMovies.push_back(i);                                        // record movie #
            userMovieRatings.push_back(users[user][i]);                     // record user rating
        }
    }

    // create current User
    User curUser(userMovieRatings, 0.0, 0, userMovieRatings);

    // find similar trainers
    for(i = 0; i < TRAINERS; i++) {                                     // for all trainers
        simMovieRatings.clear();                                        // clear movie ratings (same vector from previous trainer)
        for(j = 0; j < testNum; j++) {                                  // for all movies that were rated by the user
            if(training[i][userMovies[j]] < 1) {                        // if trainer didn't rate the movie
                simMovieRatings.push_back(0);                           // then set movie rating to 0
            } else {
                simMovieRatings.push_back(training[i][userMovies[j]]);        // get movie rating from trainer
            }
        }

        for(j = 0; j < MOVIES; j++) {
            allSimMovieRatings[j] = training[i][j];
        }

        reqRating = training[i][movie];             // get trainer's rating 
        if(reqRating < 1) {                         // if trainer's rating is < 1, set to 0
            reqRating = 0;
        }

        // create similar user
        User simUser(simMovieRatings, 0.0, reqRating, allSimMovieRatings);

        dot = 0.0;
        ru1i = 0.0;
        ru2i = 0.0;

        // find similarity
        for(j = 0; j < testNum; j++) {
            if(userMovieRatings[j] != 0 && simMovieRatings[j] != 0) {           // if movie is rated by both users
                dot += (userMovieRatings[j] * simMovieRatings[j]);              // numer dot product
                ru1i += (userMovieRatings[j] * userMovieRatings[j]);            // denom left
                ru2i += (simMovieRatings[j] * simMovieRatings[j]);              // denom right
            }
        }
        ru1i = sqrt(ru1i);
        ru2i = sqrt(ru2i);

        if(ru1i == 0 || ru2i == 0) {                            // if denom = 0
            wau = 0.0;                                          // then similarity = 0
        } else {
            wau = dot / (ru1i * ru2i);                          // else, compute similarity
        }

        simUser.setSimilarity(wau);             // set similarity b/w user and trainer

        if(reqRating != 0) {
            if(relUsersCount < KSIMILAR) {                  // if < KSIMILAR relevant trainers
                relUsers.push_back(simUser);                // record trainers ID
                relUsersCount++;                            // increment # of relevant trainers
            } else {                                        // else (already KSIMILAR trainers)
                for(k = 0; k < KSIMILAR; k++) {             // find least similar trainer
                    if(relUsers[k].getSimilarity() < relUsers[l].getSimilarity()) {
                        l = k;
                    }
                }
                relUsers[l] = simUser;                      // replace least similar trainer with current trainer
            }
        }
        
    }

    numerator = 0.0;
    denominator = 0.0;

    // user rating prediction calculation
    for(i = 0; i < relUsers.size(); i++) {                 // for k similar trainers
        wau = relUsers[i].getSimilarity();                     // get trainer's similarity
        rui = relUsers[i].getReqRating();                  // get trainer's needed rating
        numerator += (wau * (float)rui);              // calc numer
        denominator += wau;                    // calc denom
    }

    if(denominator != 0)
        rating = round(numerator / denominator);

    if(rating <= 0) {                       // if no one else rated the movie
        rating = 3;      // rating = user average
    }

    // restrict rating to a maximum
    if(rating > MAX_RATING) {
        rating = MAX_RATING;
    } 

    relUsers.clear();

    return rating;
}

int customMethod(int user, int movie) {

    int rating = 0;
    int r1 = cosineSimilarity(user, movie);
    int r2 = itemBased(user, movie);
    int r3 = pearsonCorrelation(user, movie);
    int r4 = pearsonCase(user, movie);

    rating = round((r1 + r2 + r3 + r4) / 4);

    return rating;
}

int main() {

    int i, j;
    int user;
    int movie;
    int rating;

    // open training file
    FILE *fpTrain = fopen("train.txt", "r");
    if(fpTrain == NULL) {
        fprintf(stderr, "Unable to open train.txt file\n");
        exit(1);
    }

    // put training data into ratings
    while(fscanf(fpTrain, "%d %d %d", &user, &movie, &rating) != EOF) {
        training[user-1][movie-1] = rating;
    }
    fclose(fpTrain);

    // variables for test file choice
    int fileNumber;
    int method;
    std::string fileIn;
    std::string fileOut;

    // choose test file
    printf("1 - test5.txt\n");
    printf("2 - test10.txt\n");
    printf("3 - test20.txt\n");
    scanf("%d", &fileNumber);

    switch(fileNumber) {
        case 1:
            offset = 200;
            testNum = 5;
            fileIn = "test5.txt";
            fileOut = "result5.txt";
            break;
        case 2:
            offset = 300;
            testNum = 10;
            fileIn = "test10.txt";
            fileOut = "result10.txt";
            break;
        case 3:
            offset = 400;
            testNum = 20;
            fileIn = "test20.txt";
            fileOut = "result20.txt";
            break;
    }

    // open test file
    FILE *fpTest = fopen(fileIn.c_str(), "r");
    if(fpTest == NULL) {
        fprintf(stderr, "Unable to open test file\n");
        exit(1);
    }

    // put test file data into ratings
    // 0 = needs rating 
    // 1-5 = rated
    // -1 = ignore that movie for that user
    while(fscanf(fpTest, "%d %d %d", &user, &movie, &rating) != EOF) {
        users[user-1-offset][movie-1] = rating;
    }
    fclose(fpTest);

    // calculate ratings
    printf("1 - Cosine Similarity\n");
    printf("2 - Pearson's Correlation\n");
    printf("3 - Pearson's Inverse Correlation\n");
    printf("4 - Pearson's Case Correlation\n");
    printf("5 - Item Based Collaborative Filtering\n");
    printf("6 - Custom Method\n");
    scanf("%d", &method);

    for(i = 0; i < USERS; i++) {
        for(j = 0; j < MOVIES; j++) {
            rating = -1;
            if(users[i][j] == 0) {
                switch(method){
                case 1:
                    rating = cosineSimilarity(i, j);
                    break;
                case 2:
                    rating = pearsonCorrelation(i, j);
                    break;
                case 3:
                    rating = pearsonInverse(i, j);
                    break;
                case 4:
                    rating = pearsonCase(i, j);
                    break;
                case 5:
                    rating = itemBased(i, j);
                    break;
                case 6:
                    rating = customMethod(i, j);
                    break;
                }
            }
            results[i][j] = rating;
        }
    }

    // open result file
    FILE *fpResult = fopen(fileOut.c_str(), "w");
    if(fpResult == NULL) {
        fprintf(stderr, "Unable to open result_.txt file\n");
        exit(1);
    }

    // put results into result file
    for(i = 0; i < USERS; i++) {
        for(j = 0; j < MOVIES; j++) {
            if(results[i][j] != -1) {
                user = i + offset + 1;
                movie = j + 1;
                rating = results[i][j];
                fprintf(fpResult, "%d %d %d\n", user, movie, rating);
            }
        }
    }
    fclose(fpResult);

    return 0;
}