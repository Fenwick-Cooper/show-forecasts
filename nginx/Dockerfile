FROM node:14.16.0-alpine as build

RUN apk update && apk add --update --no-cache -f git && \
    git clone --depth 1 https://github.com/mljar/mercury.git /app

WORKDIR /app/frontend

RUN yarn install --frozen-lockfile && yarn build

# The second stage
# Copy React static files and start nginx
FROM nginx:stable-alpine
COPY --from=build /app/frontend/build /usr/share/nginx/html
