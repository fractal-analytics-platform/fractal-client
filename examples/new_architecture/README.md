Follow these steps from main repo folder.


## First-time setup

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

5. Install `http` (e.g. with `sudo apt install httpie`, or with conda, snap, ...).

## Every time you want to create `data/fractal_server.db`

```
poetry run alembic revision --autogenerate -m 'init'
poetry run alembic upgrade head
```

## Every time you want to start the server and register the test user

1. `poetry run server`
2. In a different terminal: `http POST localhost:8000/auth/register email=test@me.com password=test`


Now you are ready to go, with `poetry run client` commands.

## Or you can use a few scripts..

```bash
./server_setup.sh
./server_start.sh
# Open new terminal
./register_user.sh
poetry run client project list
```
