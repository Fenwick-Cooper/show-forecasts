#!/bin/sh

until cd ${APP_DIR}; do
    echo "Waiting for server volume..."
done

until . ${VENV_DIR}/bin/activate; do
    echo "activating virtual environment..."
done

until cd ${WORK_HOME}/mercury &&
    python manage.py migrate; do
    echo "Waiting for db to be ready..."
    sleep 2
done

echo "Docker: Add mercury web applications"
# for filepath in ${APP_DIR}/webapps/*.ipynb; do
#     echo "Docker: Add " $filepath
#     python manage.py add $filepath
# done
echo "Docker: Add " ${APP_DIR}/webapps/ecmwf-ensemble-visualization.ipynb
python manage.py add ${APP_DIR}/webapps/ecmwf-ensemble-visualization.ipynb

echo "Docker: Collect statics, for Admin Panel, DRF views"
python manage.py collectstatic --noinput

echo "Docker: Try to create super user, if doesnt exist"
python manage.py createsuperuser --noinput

echo "Docker: Start worker and beat service"
celery -A server worker --loglevel=info -P gevent --concurrency 4 -E -Q celery,ws &
celery -A server beat --loglevel=error --max-interval 60 &

echo "Docker: Start daphne server"
daphne server.asgi:application --bind 0.0.0.0 --port 9000

#gunicorn server.wsgi --bind 0.0.0.0:8000 --workers 4 --threads 4

# for debug
#python manage.py runserver 0.0.0.0:9000
