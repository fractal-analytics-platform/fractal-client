Follow these steps from main repo folder.


## THE FIRST TIME ONLY

1. The first time, run
```
git fetch
git checkout client
git pull
poetry install
```

2. Write the following content in `.fractal_server.env`:
```
DEPLOYMENT_TYPE=testing
JWT_SECRET_KEY=secret
DATA_DIR_ROOT=/tmp/
SQLITE_PATH=./data/fractal_server.db
```

3. Write the following content in `.fractal.env`:
```
FRACTAL_USER=test@me.com
FRACTAL_PASSWORD=test
```

4.
```
mkdir data
mkdir fractal/server/migrations/versions/
```

5. `sudo apt install httpie`

# EVERY TIME YOU NEED TO CREATE THE DB

```
poetry run alembic revision --autogenerate -m 'init'
poetry run alembic upgrade head
```

# TO START THE SERVER

1. `poetry run server`
2. In a different terminal: `http POST localhost:8000/auth/register email=test@me.com password=test`


Ready to go :)
