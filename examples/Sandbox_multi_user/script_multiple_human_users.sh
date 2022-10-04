# Before anything else happens, someone (AKA Joel/Gustavo) should register the
# user test@me.com. After this command, any other repetition of the same yields
# "REGISTER_USER_ALREADY_EXISTS" message and does nothing.
curl -d '{"email":"test@me.com", "password":"test"}' -H "Content-Type: application/json" -X POST localhost:8000/auth/register
echo

# Let's create a folder for each human user, including a .fractal.env file
# (this is the simplest way, so that the user does not need to enter
# user/password at each command)
rm -r human-user-1 human-user-2
mkdir human-user-1
mkdir human-user-2
echo -e "FRACTAL_USER=test@me.com\nFRACTAL_PASSWORD=test" > human-user-1/.fractal.env
echo -e "FRACTAL_USER=test@me.com\nFRACTAL_PASSWORD=test" > human-user-2/.fractal.env

# This block is what human user 1 would do
cd human-user-1                          # Enter your own folder
export FRACTAL_CACHE_PATH=cache          # Set the cache folder (notice that this points to human-user-1/cache)
fractal project list > /dev/null 2>&1    # Run some fractal client command
echo
cd ..

# This block is what human user 2 would do
cd human-user-2                          # Enter your own folder
export FRACTAL_CACHE_PATH=cache          # Set the cache folder (notice that this points to human-user-2/cache)
fractal project list > /dev/null 2>&1    # Run some fractal client command
echo
cd ..

# We can check that both user have a token available in their own cache folder,
# which is used to authenticate their fractal client commands. The tokens can be
# the same or different, because each token has a finite lifetime and then is
# replaced by a new one
for F in `ls human-user-*/cache/session`; do
    echo $F
    cat $F
    echo
    echo
done
